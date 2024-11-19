#@Time      2024/10/30 16:30
#@Author    hzb
#将fq文件的提取kmer与计算矩阵、计算相似度、绘出发育树进行整合
#11.5号更改，添加对于原始文件的更改，对于原文件不进行过滤频率

import subprocess
import os
import glob
import numpy as np
import pandas as pd
from scipy import stats

#对进行计算后的头部kmer数据进行相似度分析

# 主程序

# 2. 计算Jaccard指数
class similary:#q1:mum3 have 2 hashvalues
    def __init__(self,list1,list2,k):
        #print("step1=====",type(list1),list1)
        self.k = k
        self.minhash1 = self.getkmer(list1)
        self.minhash2= self.getkmer(list2)
        self.fre1,self.fre2 = self.getfrequence(list1,list2)
        #print(fre1,fre2)
        #print("step2====",self.minhash1)
    def getkmer(self,list_dict):
        strings = []
        for t in list_dict:
            strings.append(t[0])
        return strings

    def getfrequence(self,set1,set2):
        # 获取所有 k-mer 的集合
        # 转换为字典以便于查找
        dict1 = {kmer: int(freq) for kmer, freq in set1}
        dict2 = {kmer: int(freq) for kmer, freq in set2}
        # 合并频率数列
        merged_frequencies1 = []
        merged_frequencies2 = []
        # 计算频率
        for kmer in dict1.keys() | dict2.keys():  # 使用集合的并集
            freq1 = dict1.get(kmer, 0)  # 从 dict1 获取频率，没有则为 0
            freq2 = dict2.get(kmer, 0)  # 从 dict2 获取频率，没有则为 0
            merged_frequencies1.append(freq1)
            merged_frequencies2.append(freq2)
        return merged_frequencies1,merged_frequencies2

    def mash_distance(self) -> float:
        jaccard = self.jaccard_distance()
        if jaccard == 0:
            mash_distance = 0
        else:
            mash_distance = -1 / self.k * np.log(2 * jaccard / (1 + jaccard))
        #print("mashdistance", mash_distance)
        return  mash_distance
    def jaccard_distance(self):
        #print('------',self.minhash1,"\n",self.minhash2)

        a = set(self.minhash1)
        b = set(self.minhash2)
        #print('step3===',a,b)
        intersection = a.intersection(b)
        union = a.union(b)
        c = len(intersection)
        d = len(union)
        #sim = fr'{c}/{d}'#用来得到匹配数量
        #print('he')
        #calculate and obtain the similarity
        jaccard_index = len(intersection) / len(union) if len(union) != 0 else 0

        return jaccard_index

    # pearson corrlation
    def pearson(self):
        correlation, p_value = stats.pearsonr(self.fre1, self.fre2)
        #orr, p = stats.spearmanr(self.freq_list1, self.freq_list2)
        # 输出结果
        #print("Pearson 相关系数:", correlation)
        #print("p-value:", p_value)
        # print("spearman 相关系数:", correlation)
        # print("p-value:", p_value)
        return correlation

    # 余弦相似度
    def cosine_similarity(self):
        # 计算向量的点积
        dot_product = np.dot(self.fre1, self.fre2)
        # 计算向量的范数（即向量的长度）
        norm_vector1 = np.linalg.norm(self.fre1)
        norm_vector2 = np.linalg.norm(self.fre2)
        # 计算余弦相似性
        cosine_sim = round(dot_product / (norm_vector1 * norm_vector2),3)
        #print('余弦相似度:', cosine_sim)
        return cosine_sim

    # 欧几里得距离
    def euclidean_distance(self):
        a = np.array(self.fre1)
        b = np.array(self.fre2)
        c = np.linalg.norm(a - b)
        #print("欧几里得距离:", c)
        return c


def get_kmer(input_dir,output_file,top_n,kmer_size,cutoff,root_node):
    # 设置文件路径和参数

    #kmer_size = 16
    threads = 5
    #top_n = 5000
    #cutoff = 2

    # 查找所有 .fq 文件
    fq_files1 = glob.glob(os.path.join(input_dir, "*.fq"))
    fq_files2 = glob.glob(os.path.join(input_dir, "*.fasta"))
    fq_files = fq_files1+fq_files2
    matrix = []
    row_file = ''
    # 处理每个 .fq 文件
    for file in fq_files:
        print(f"Processing {file}...")

        # 获取文件名，不带路径和扩展名，如果文件名等于原始文件名，就不过滤频率
        row_name = os.path.basename(root_node).replace("*.{fq,fasta}", "")
        filename = os.path.basename(file).replace("*.{fq,fasta}", "")
        print("----------------",row_name,filename)
        # 定义输出文件名
        output_top = os.path.join(output_file, f"{filename}_top_kmers.txt")
        print("----------------",output_top)
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
            result = None
            result = subprocess.run(command, stdout=outfile, stderr=subprocess.PIPE, text=True)
    # # 打印标准错误输出（如果存在）
    #if result.stderr:
         #print("错误输出:", result.stderr)
    return row_file

# 1. 数据提取
def extract_kmers(file_path, num_lines):
    """
    从指定文件中提取最后num_lines行的k-mer及其频率。
    参数：
        file_path (str): 文件路径。
        num_lines (int): 提取的行数，默认为1000。
    返回：
        set: 包含k-mer的集合。
    """
    #print(file_path)
    with open(file_path, 'r') as file:
        lines = file.readlines()# 读取所有行
        last_lines = lines[-num_lines:]# 提取最后num_lines行
        #print(last_lines)
        kmers = set()
        for line in last_lines:# 提取k-mer（假设每行以tab分隔，第一列为k-mer）
            parts = line.strip().split('\t')
            if len(parts) >= 1:
                kmers.add((parts[0],parts[2]))
    #print(kmers)
    return kmers

def load_all_kmers(directory, num_lines):
    """
       加载目录下所有txt文件的k-mer集合。
       参数：
           directory (str): 目录路径。
           num_lines (int): 提取的行数，默认为1000。
       返回：
           dict: 键为文件名，值为k-mer集合。
       """
    kmers_dict = {}
    #print("ddd____",directory)
    for filename in os.listdir(directory):
        if filename.endswith('.txt'):
            file_path = os.path.join(directory, filename)#get file path
            kmers = extract_kmers(file_path, num_lines)##get kmer
            kmers_dict[filename] = kmers#save to dict
    return kmers_dict


# 3. 构建距离矩阵
def build_distance_matrix(kmers_dict,kmer_size):
    """
        构建基于Jaccard距离的距离矩阵。
        参数：
            kmers_dict (dict): 键为文件名，值为k-mer集合。
        返回：
            pd.DataFrame: 距离矩阵。
        """
    files = list(kmers_dict.keys())
    num_files = len(files)#obtainthe number of included files
    #print(num_files)
    distance_matrix = pd.DataFrame(np.zeros((num_files, num_files)), index=files, columns=files)
    #print('distance_matrix',distance_matrix)
    for i, file1 in enumerate(files):
        for j, file2 in enumerate(files):
            if i < j:
                #print(file1)
                set1 = kmers_dict[file1]
                set2 = kmers_dict[file2]
                sim = similary(set1,set2,kmer_size)
                jaccard = sim.jaccard_distance()
                #print(i,j,'jaccard',jaccard)
                mash = sim.mash_distance()
                r = sim.pearson()
                cos = sim.cosine_similarity()
                ou = sim.euclidean_distance()
                #print("pearson",r)
                #distance =fr'{jaccard}cos:{cos}'#(1-jaccard)#(mash*(1-cos))**0.5
                distance = mash#1-cos#(mash*(1-cos))**0.5

                #fill the matrix
                distance_matrix.at[file1, file2] = distance
                distance_matrix.at[file2, file1] = distance
            elif i == j:
                distance_matrix.at[file1, file2] = 0.0

    # 创建 DataFrame
    df = pd.DataFrame(distance_matrix)
    # 导出到 CSV 文件
    df.to_csv('distance_matrix-mashcos-seed42.csv', index=True, header=True)
    return distance_matrix

def similary_distance(directory,top_kmer,kmer_size):
    # 设置数据目录
    #directory = '/mnt/00.zhibo/6M100Retrs/bing/mutated2'  # 请替换为您的数据目录路径
    # 数据提取
    kmers_data = load_all_kmers(directory, top_kmer)
    print("已提取所有文件的k-mer数据。")
    # 构建距离矩阵
    distance_df = build_distance_matrix(kmers_data,kmer_size)
    print("已构建Jaccard距离矩阵。")
    print(distance_df)
    return distance_df


def main(input_dir = "/mnt/00.zhibo/random/random_t1/",output_file =  "/mnt/00.zhibo/random/random_t1"):
    #input_dir = "/mnt/00.zhibo/random/row_min_t/"#该地址为原始fq文件集合的文件夹/mnt/path/
    #output_file = "/mnt/00.zhibo/random/row_min"#default="/mnt/path2"
    print(output_file)
    top_kmer = 5000
    kmer_size = 31
    cutoff = '2'#过滤出现频率
    root_node = 'row_min_.fasta'#是否设置发育树根节点#default:name.txt
    #计算提取fq文件的头部kmer
    root_node2 = get_kmer(input_dir,output_file,top_kmer,kmer_size,cutoff,root_node)#输入文件输出文件抽样比例k值过滤频率根节点
    distance = similary_distance(output_file,top_kmer,kmer_size)
    # root_node = 'row_random_5000mmh3_kmers.txt'
    return distance,root_node2

if __name__ == '__main__':
    main()