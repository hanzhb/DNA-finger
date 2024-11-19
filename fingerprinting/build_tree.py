# 4. 构建系统发育树
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from matplotlib import pyplot as plt


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
    plt.savefig('kc.png')

