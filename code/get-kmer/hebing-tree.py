#@Time      2024/10/30 16:30
#@Author    hzb
#将fq文件的提取kmer与计算矩阵、计算相似度、绘出发育树进行整合
#11.5号更改，添加对于原始文件的更改，对于原文件不进行过滤频率
import subprocess
import os
import glob
import tree_1023

def get_kmer(input_dir,output_file,top_n,kmer_size,cutoff,root_node):
    # 设置文件路径和参数

    #kmer_size = 16
    threads = 5
    #top_n = 5000
    #cutoff = 2

    # 查找所有 .fq 文件
    fq_files = glob.glob(os.path.join(input_dir, "*"))
    matrix = []
    row_file = ''
    # 处理每个 .fq 文件
    for file in fq_files:
        print(f"Processing {file}...")

        # 获取文件名，不带路径和扩展名，如果文件名等于原始文件名，就不过滤频率
        row_name = os.path.basename(root_node).replace(".*", "")
        filename = os.path.basename(file).replace(".*", "")
        # 定义输出文件名
        output_top = os.path.join(output_file, f"{filename}_5000mmh3_kmers.txt")

        if filename == row_name:
            row_file = os.path.basename(output_top)
            # 构建 kc-c4-100 命令
            command = [
                "/mnt/00.zhibo/6M100Retrs/bing/kc-c4-100",
                "-k", str(kmer_size),
                "-b", "1000000",
                "-t", str(threads),
                "-N", str(top_n),
                "-c", "",
                "-o", "/mnt/00.zhibo/6M100Retrs/bing/test2/top1.txt",
                file,
                ">",  # 用于重定向输出的符号，Python subprocess 不支持
                output_top
            ]
        else:
            # 构建 kc-c4-100 命令
            command = [
                "/mnt/00.zhibo/6M100Retrs/bing/kc-c4-100",
                "-k", str(kmer_size),
                "-b", "1000000",
                "-t", str(threads),
                "-N", str(top_n),
                "-c", str(cutoff),
                "-o", "/mnt/00.zhibo/6M100Retrs/bing/test2/top1.txt",
                file,
                ">",  # 用于重定向输出的符号，Python subprocess 不支持
                output_top
            ]
        # 因为 Python 的 subprocess 不支持 shell 重定向，采用直接输出到文件的方式
        with open(output_top, 'w') as outfile:
        # 调用C程序并捕获其标准输出
            result = subprocess.run(command, stdout=outfile, stderr=subprocess.PIPE, text=True)
    # 打印标准错误输出（如果存在）
    if result.stderr:
        print("错误输出:", result.stderr)
    return row_file

def main():
    input_dir = "/mnt/00.zhibo/random/row_min_t/"#该地址为原始fq文件集合的文件夹/mnt/path/
    output_file = "/mnt/00.zhibo/random/row_min"#default="/mnt/path2"
    top_kmer = 5000
    kmer_size = 31
    cutoff = '2'#过滤出现频率
    root_node = 'row_min.fasta'#是否设置发育树根节点#default:name.txt
    #计算提取fq文件的头部kmer
    root_node = get_kmer(input_dir,output_file,top_kmer,kmer_size,cutoff,root_node)#输入文件输出文件抽样比例k值过滤频率根节点
    #调用tree方法，获得系统发育树
    #root_node = 'row_random_5000mmh3_kmers.txt'
    tree_1023.main(output_file,top_kmer,kmer_size,root_node)

if __name__ == '__main__':
    main()



