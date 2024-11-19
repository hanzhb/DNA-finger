#@Time      2024/11/11 16:30
#@Author    hzb
#

import  build_tree
import similary_distance
import mutated_copy
import gunzipfile

input_row_fasta = '/mnt/00.zhibo/random/row_min.fasta'
output_file= "/mnt/00.zhibo/random/random_t3/"

mutated_copy.bingxing(input_row_fasta,output_file,3)
#创建随机文件fq，循环几次代表几级串行拷贝
mutated_copy.chuanxing(input_row_fasta,output_file,cycle=6)#5 times copy
#
# #解压文件
gunzipfile.decompress_and_move_gz_files(output_file,output_file)

#计算距离矩阵
distance,nood = similary_distance.main(output_file,output_file)
#build tree
build_tree.build_phylogenetic_tree(distance,nood)