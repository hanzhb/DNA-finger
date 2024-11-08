#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>
#include <stdint.h>
#include <stdbool.h>
#include <ctype.h>
#include <zlib.h>
#include <time.h>
#include <unistd.h>

#define MAX_KMER_SIZE 31   // Maximum k-mer length
#define HASH_TABLE_SIZE 10000019  // A large prime number for hash table size
#define BUFFER_SIZE 4096    // Buffer size for file reading

// Structure to hold k-mer counts
typedef struct {
    uint64_t kmer;
    uint32_t count;
} KmerCount;

// Node for handling collisions in hash table
typedef struct KmerNode {
    uint64_t kmer;
    uint32_t count;
    struct KmerNode *next;
} KmerNode;

// Thread data structure
typedef struct {
    char *sequence;
    size_t start;
    size_t end;
    int k;
    KmerNode **hash_table;
    pthread_mutex_t *mutexes;
} ThreadData;

// MurmurHash3 64-bit hash function
uint64_t murmur3_64(uint64_t key) {
    key ^= key >> 33;
    key *= 0xff51afd7ed558ccdULL;
    key ^= key >> 33;
    key *= 0xc4ceb9fe1a85ec53ULL;
    key ^= key >> 33;
    return key;
}

// Thomas Wang's 64-bit hash function
uint64_t wang_hash64(uint64_t key) {
    key = (~key) + (key << 21); // key = (key << 21) - key - 1;
    key = key ^ (key >> 24);
    key = (key + (key << 3)) + (key << 8); // key * 265
    key = key ^ (key >> 14);
    key = (key + (key << 2)) + (key << 4); // key * 21
    key = key ^ (key >> 28);
    key = key + (key << 31);
    return key;
}

// Function pointer for the chosen hash function
uint64_t (*hash_function)(uint64_t);

// Convert DNA sequence to 2-bit representation
uint64_t dna_to_int(const char *seq, int len) {
    uint64_t result = 0;
    for (int i = 0; i < len; i++) {
        result <<= 2;
        switch (seq[i]) {
            case 'A': case 'a': result |= 0; break;
            case 'C': case 'c': result |= 1; break;
            case 'G': case 'g': result |= 2; break;
            case 'T': case 't': result |= 3; break;
            default: return UINT64_MAX; // Non-ACGT character
        }
    }
    return result;
}

// Initialize hash table
KmerNode **init_hash_table() {
    KmerNode **table = calloc(HASH_TABLE_SIZE, sizeof(KmerNode *));
    if (!table) {
        perror("Hash table allocation failed");
        exit(EXIT_FAILURE);
    }
    return table;
}

// Insert or increment k-mer count
void insert_kmer(KmerNode **table, pthread_mutex_t *mutexes, uint64_t kmer) {
    uint64_t hash = hash_function(kmer) % HASH_TABLE_SIZE;
    pthread_mutex_lock(&mutexes[hash % 256]); // Use a subset of mutexes
    KmerNode *node = table[hash];
    while (node) {
        if (node->kmer == kmer) {
            node->count++;
            pthread_mutex_unlock(&mutexes[hash % 256]);
            return;
        }
        node = node->next;
    }
    // K-mer not found, add new node
    KmerNode *new_node = malloc(sizeof(KmerNode));
    if (!new_node) {
        perror("KmerNode allocation failed");
        pthread_mutex_unlock(&mutexes[hash % 256]);
        exit(EXIT_FAILURE);
    }
    new_node->kmer = kmer;
    new_node->count = 1;
    new_node->next = table[hash];
    table[hash] = new_node;
    pthread_mutex_unlock(&mutexes[hash % 256]);
}

// Determine file format: FASTA or FASTQ
typedef enum { UNKNOWN, FASTA, FASTQ } FileFormat;

FileFormat get_file_format(const char *filename) {
    const char *dot = strrchr(filename, '.');
    if (dot) {
        if (strcmp(dot, ".fa") == 0 || strcmp(dot, ".fasta") == 0)
            return FASTA;
        if (strcmp(dot, ".fq") == 0 || strcmp(dot, ".fastq") == 0)
            return FASTQ;
        if (strcmp(dot, ".gz") == 0) {
            // Check the extension before .gz
            size_t len = strlen(filename);
            for (size_t i = len - 4; i > 0; i--) {
                if (filename[i] == '.') {
                    if (strncmp(&filename[i], ".fa", 3) == 0 || strncmp(&filename[i], ".fasta", 6) == 0)
                        return FASTA;
                    if (strncmp(&filename[i], ".fq", 3) == 0 || strncmp(&filename[i], ".fastq", 6) == 0)
                        return FASTQ;
                    break;
                }
            }
        }
    }
    return UNKNOWN;
}

// Read sequences from file into a single string
char *read_sequences(const char *filename, size_t *total_length) {
    gzFile file = gzopen(filename, "r");
    if (!file) {
        perror("Error opening file");
        return NULL;
    }

    FileFormat format = get_file_format(filename);
    if (format == UNKNOWN) {
        fprintf(stderr, "Unknown file format\n");
        gzclose(file);
        return NULL;
    }

    char *sequences = NULL;
    size_t capacity = BUFFER_SIZE;
    size_t length = 0;
    char buffer[BUFFER_SIZE];
    int bytes_read;
    bool in_sequence = false;
    int line_part = 0; // For FASTQ format

    sequences = malloc(capacity);
    if (!sequences) {
        perror("Memory allocation error");
        gzclose(file);
        return NULL;
    }

    while ((bytes_read = gzread(file, buffer, BUFFER_SIZE)) > 0) {
        for (int i = 0; i < bytes_read; i++) {
            char c = buffer[i];
            if (format == FASTA) {
                if (c == '>') {
                    in_sequence = false;
                    // Skip the header line
                    while (i < bytes_read && buffer[i] != '\n') i++;
                } else if (c == '\n') {
                    in_sequence = true;
                } else if (in_sequence && isalpha(c)) {
                    if (length + 1 >= capacity) {
                        capacity *= 2;
                        sequences = realloc(sequences, capacity);
                        if (!sequences) {
                            perror("Memory reallocation error");
                            gzclose(file);
                            return NULL;
                        }
                    }
                    sequences[length++] = toupper(c);
                }
            } else if (format == FASTQ) {
                // FASTQ format: sequence lines are every 4th line starting from the second line
                if (c == '\n') {
                    line_part = (line_part + 1) % 4;
                } else if (line_part == 1 && isalpha(c)) {
                    if (length + 1 >= capacity) {
                        capacity *= 2;
                        sequences = realloc(sequences, capacity);
                        if (!sequences) {
                            perror("Memory reallocation error");
                            gzclose(file);
                            return NULL;
                        }
                    }
                    sequences[length++] = toupper(c);
                }
            }
        }
    }

    if (bytes_read < 0) {
        fprintf(stderr, "Error reading file\n");
        free(sequences);
        gzclose(file);
        return NULL;
    }

    sequences[length] = '\0';
    *total_length = length;
    gzclose(file);
    return sequences;
}

void *count_kmers(void *arg) {
    ThreadData *data = (ThreadData *)arg;
    const char *seq = data->sequence;
    int k = data->k;
    size_t start = data->start;
    size_t end = data->end;
    KmerNode **hash_table = data->hash_table;
    pthread_mutex_t *mutexes = data->mutexes;

    size_t processed = 0;
    for (size_t i = start; i <= end - k; i++) {
        uint64_t kmer = dna_to_int(&seq[i], k);
        if (kmer == UINT64_MAX) continue; // Skip invalid k-mers
        insert_kmer(hash_table, mutexes, kmer);
        processed++;
        // Optional: progress reporting
        if (processed % 1000000 == 0) {
            printf("Thread %ld processed %zu k-mers\n", pthread_self(), processed);
        }
    }
    return NULL;
}

// Comparator function for qsort (sort by hash value in ascending order)
int compare_kmer_hashes(const void *a, const void *b) {
    const KmerCount *kmer_a = (const KmerCount *)a;
    const KmerCount *kmer_b = (const KmerCount *)b;
    uint64_t hash_a = hash_function(kmer_a->kmer);
    uint64_t hash_b = hash_function(kmer_b->kmer);
    if (hash_a < hash_b) return -1;
    else if (hash_a > hash_b) return 1;
    else return 0;
}

void print_usage(const char *program_name) {
    fprintf(stderr, "Usage: %s [OPTIONS] <sequence_file>\n", program_name);
    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  -k <int>        K-mer size (required)\n");
    fprintf(stderr, "  -t <int>        Number of threads (default: 1)\n");
    fprintf(stderr, "  -o <file>       Output file for all k-mer counts\n");
    fprintf(stderr, "  -F <file>       Output file for top N k-mers with smallest hash values\n");
    fprintf(stderr, "  -N <int>        Number of top k-mers with smallest hash values to output (default: 1000)\n");
    fprintf(stderr, "  -c <int>        Coverage cutoff for k-mers (default: 0)\n");
    fprintf(stderr, "  -w              Use Wang's hash function (default: MurmurHash3)\n");
    fprintf(stderr, "  -h              Display this help message and exit\n");
}

int main(int argc, char *argv[]) {
    int k = 0;
    int num_threads = 1;
    const char *sequence_file = NULL;
    const char *output_file = NULL;
    const char *F_output_file = NULL; // Output file for top N k-mers
    int N = 1000; // Default number of k-mers to output
    int coverage_cutoff = 0; // Coverage cutoff
    bool use_wang_hash = false; // Flag to use Wang's hash function
    int opt;

    // Parse command-line options
    while ((opt = getopt(argc, argv, "k:t:o:F:N:c:wh")) != -1) {
        switch (opt) {
            case 'k':
                k = atoi(optarg);
                break;
            case 't':
                num_threads = atoi(optarg);
                break;
            case 'o':
                output_file = optarg;
                break;
            case 'F':
                F_output_file = optarg;
                break;
            case 'N':
                N = atoi(optarg);
                break;
            case 'c':
                coverage_cutoff = atoi(optarg);
                break;
            case 'w':
                use_wang_hash = true;
                break;
            case 'h':
                print_usage(argv[0]);
                return EXIT_SUCCESS;
            default:
                print_usage(argv[0]);
                return EXIT_FAILURE;
        }
    }

    // Check if required arguments are provided
    if (optind >= argc || k <= 0 || k > MAX_KMER_SIZE || num_threads <= 0) {
        print_usage(argv[0]);
        return EXIT_FAILURE;
    }

    sequence_file = argv[optind];

    if (k <= 0 || k > MAX_KMER_SIZE) {
        fprintf(stderr, "Invalid k-mer size. Must be between 1 and %d.\n", MAX_KMER_SIZE);
        return EXIT_FAILURE;
    }
    if (num_threads <= 0) {
        fprintf(stderr, "Invalid number of threads.\n");
        return EXIT_FAILURE;
    }

    // Set the hash function based on the user's choice
    hash_function = use_wang_hash ? wang_hash64 : murmur3_64;

    // Read sequences from file
    size_t total_length;
    char *sequence = read_sequences(sequence_file, &total_length);
    if (!sequence) {
        fprintf(stderr, "Failed to read sequences.\n");
        return EXIT_FAILURE;
    }
    printf("Total sequence length: %zu\n", total_length);

    // Initialize hash table and mutexes
    KmerNode **hash_table = init_hash_table();
    pthread_mutex_t mutexes[256]; // Use 256 mutexes for simplicity
    for (int i = 0; i < 256; i++) {
        pthread_mutex_init(&mutexes[i], NULL);
    }

    // Start timing
    clock_t start = clock();

    // Create threads
    pthread_t threads[num_threads];
    ThreadData thread_data[num_threads];
    size_t chunk_size = total_length / num_threads;

    for (int i = 0; i < num_threads; i++) {
        thread_data[i].sequence = sequence;
        thread_data[i].k = k;
        thread_data[i].hash_table = hash_table;
        thread_data[i].mutexes = mutexes;
        thread_data[i].start = i * chunk_size;
        thread_data[i].end = (i == num_threads - 1) ? total_length - 1 : (i + 1) * chunk_size + k - 1;

        if (pthread_create(&threads[i], NULL, count_kmers, &thread_data[i]) != 0) {
            perror("Thread creation failed");
            return EXIT_FAILURE;
        }
    }

    // Wait for threads to finish
    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    // Stop timing
    clock_t end = clock();
    double cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC;

    printf("K-mer counting completed in %.2f seconds\n", cpu_time_used);

    // Output results to file if specified
    if (output_file) {
        FILE *out = fopen(output_file, "w");
        if (!out) {
            perror("Error opening output file");
            return EXIT_FAILURE;
        }
        for (size_t i = 0; i < HASH_TABLE_SIZE; i++) {
            KmerNode *node = hash_table[i];
            while (node) {
                // Convert integer k-mer back to DNA sequence
                char kmer_seq[MAX_KMER_SIZE + 1];
                uint64_t temp = node->kmer;
                for (int j = k - 1; j >= 0; j--) {
                    int base = temp & 0x3;
                    kmer_seq[j] = "ACGT"[base];
                    temp >>= 2;
                }
                kmer_seq[k] = '\0';
                fprintf(out, "%s\t%u\n", kmer_seq, node->count);
                node = node->next;
            }
        }
        fclose(out);
        printf("K-mer counts written to %s\n", output_file);
    } else {
        printf("No output file specified. K-mer counts not written to file.\n");
    }

    // If -F option is provided, output top N k-mers with smallest hash values
    if (F_output_file) {
        KmerCount *top_kmers = malloc(N * sizeof(KmerCount));
        if (!top_kmers) {
            perror("Memory allocation failed for top_kmers");
            return EXIT_FAILURE;
        }
        int top_count = 0;

        for (size_t i = 0; i < HASH_TABLE_SIZE; i++) {
            KmerNode *node = hash_table[i];
            while (node) {
                if (node->count > coverage_cutoff) {
                    uint64_t hash = hash_function(node->kmer);
                    if (top_count < N) {
                        top_kmers[top_count].kmer = node->kmer;
                        top_kmers[top_count].count = node->count;
                        top_count++;
                        if (top_count == N) {
                            qsort(top_kmers, N, sizeof(KmerCount), compare_kmer_hashes);
                        }
                    } else if (hash < hash_function(top_kmers[N-1].kmer)) {
                        top_kmers[N-1].kmer = node->kmer;
                        top_kmers[N-1].count = node->count;
                        qsort(top_kmers, N, sizeof(KmerCount), compare_kmer_hashes);
                    }
                }
                node = node->next;
            }
        }

        // Output the top N k-mers to the specified file
        FILE *F_out = fopen(F_output_file, "w");
        if (!F_out) {
            perror("Error opening output file for top N k-mers with smallest hash values");
            free(top_kmers);
            return EXIT_FAILURE;
        }

        for (int i = 0; i < top_count; i++) {
            // Convert integer k-mer back to DNA sequence
            char kmer_seq[MAX_KMER_SIZE + 1];
            uint64_t temp = top_kmers[i].kmer;
            for (int j = k - 1; j >= 0; j--) {
                int base = temp & 0x3;
                kmer_seq[j] = "ACGT"[base];
                temp >>= 2;
            }
            kmer_seq[k] = '\0';
            uint64_t hash = hash_function(top_kmers[i].kmer);
            fprintf(F_out, "%s\t%u\t%llu\n", kmer_seq, top_kmers[i].count, (unsigned long long)hash);
        }
        fclose(F_out);
        printf("Top %d k-mers with smallest hash values written to %s using %s hash function\n", 
               top_count, F_output_file, use_wang_hash ? "Wang's" : "MurmurHash3");

        free(top_kmers);
    }

    // Clean up
    free(sequence);
    for (size_t i = 0; i < HASH_TABLE_SIZE; i++) {
        KmerNode *node = hash_table[i];
        while (node) {
            KmerNode *tmp = node;
            node = node->next;
            free(tmp);
        }
    }
    free(hash_table);
    for (int i = 0; i < 256; i++) {
        pthread_mutex_destroy(&mutexes[i]);
    }

    return EXIT_SUCCESS;
}