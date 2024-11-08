#对距离结果进行正常建树，加上根节点，--->更改cos的计算方法
# @Author  ASH
#@Time  2024.10.21

import csv
import os
from io import StringIO
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
import numpy as np
import pandas as pd
from itertools import combinations
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from Bio.Phylo.TreeConstruction import _DistanceMatrix
from Bio import Phylo
from scipy import stats
from sklearn.metrics import jaccard_score
from scipy.cluster.hierarchy import linkage, dendrogram, to_tree
from scipy.spatial.distance import squareform
import matplotlib.pyplot as plt
#对进行计算后的头部kmer数据进行相似度分析，构建系统进化树

# 主程序
def main(directory,top_kmer,kmer_size,root_node):
    # 设置数据目录
    #directory = '/mnt/00.zhibo/6M100Retrs/bing/mutated2'  # 请替换为您的数据目录路径

    # 数据提取
    kmers_data = load_all_kmers(directory, top_kmer)
    print("已提取所有文件的k-mer数据。")

    # 构建距离矩阵
    distance_df = build_distance_matrix(kmers_data,kmer_size)
    print("已构建Jaccard距离矩阵。")
    print(distance_df)
    # 构建系统发育树
    linkage_matrix = build_phylogenetic_tree(distance_df,root_node)
    print("已构建系统发育树的链接矩阵。")

    # 可视化系统发育树
    #plot_dendrogram(linkage_matrix, labels=distance_df.index.tolist(),
                    #title='Phylogenetic Tree Based on Jaccard Distance')


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


# 2. 计算Jaccard指数
class similary:#q1:mum3 have 2 hashvalues
    def __init__(self,list1,list2,k):
        #print("step1=====",type(list1),list1)
        self.k = k
        self.minhash1 = self.getkmer(list1)
        self.minhash2= self.getkmer(list2)
        self.fre1,self.fre2 = self.getfre(list1,list2)
        #print(fre1,fre2)
        #print("step2====",self.minhash1)
    def getkmer(self,list_dict):
        strings = []
        for t in list_dict:
            strings.append(t[0])
        return strings

    def getfre(self,set1,set2):
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
                distance = 1-cos#(mash*(1-cos))**0.5

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


# 4. 构建系统发育树
def build_phylogenetic_tree(df_distance,root_node=''):
    labels = df_distance.index.tolist()
    df = df_distance

    # 3. 提取下三角矩阵
    # BioPython 的 DistanceMatrix 要求矩阵为下三角形式
    distance_matrix = []
    for i in range(len(df)):
        row = df.iloc[i, :i + 1].tolist()  # 提取每行的下三角部分
        distance_matrix.append(row)
    # 3. 构建BioPython的 DistanceMatrix 对象
    dm = DistanceMatrix(names=labels, matrix=distance_matrix)

    # 4. 使用UPGMA算法构建系统发育树 (你也可以选择 'nj' 近邻算法)
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)  # 也可以使用 constructor.nj(dm) 来选择近邻法

    # 5. 为树设置根节点 (如果需要指定根，可以使用下面的命令设置一个标签作为根)
    if root_node :  # 检查是否提供了根节点
        tree.root_with_outgroup({'name': root_node})  # 替换 'Your_Root_Label' 为你想指定的根节点

    # 7. 绘制系统发育树，去除内部节点标签并调整字体
    def label_func(clade):
        if clade.is_terminal():  # 仅显示叶节点（终端节点）的标签
            return clade.name
        return None  # 不显示内部节点

    # 6. 绘制有根树
    plt.figure(figsize=(12, 8))

    Phylo.draw(tree, do_show=False,  # 避免立即显示，先应用自定义样式
               label_func=label_func,  # 去除内部节点标签
               axes=plt.gca(),  # 将绘图绑定到当前图表
               show_confidence=False)  # 隐藏置信度值 # 使用 BioPython 绘制系统发育树
    plt.savefig('mutated.png')


if __name__ == "__main__":
    main()


