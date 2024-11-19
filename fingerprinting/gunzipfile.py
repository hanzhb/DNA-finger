import os
import subprocess
import shutil

def decompress_and_move_gz_files(root_dir, dest_dir):
    # 检查目标目录是否存在，如果不存在则创建
    os.makedirs(dest_dir, exist_ok=True)

    # 遍历 root_dir 下的所有文件和子文件夹
    for dirpath, _, filenames in os.walk(root_dir):
        folder_name = os.path.basename(dirpath)  # 获取当前文件夹名称
        for filename in filenames:
            if filename.endswith('.gz'):
                gz_file_path = os.path.join(dirpath, filename)
                target_filename = f"{folder_name}_{filename.replace('.gz', '.fq')}"  # 使用唯一的解压文件名
                temp_target_path = os.path.join(dirpath, target_filename)

                # 解压 .gz 文件
                try:
                    with open(temp_target_path, 'wb') as f_out:
                        subprocess.run(['gunzip', '-c', gz_file_path], check=True, stdout=f_out)
                    #os.remove(gz_file_path)  # 删除原始 .gz 文件
                    # 移动解压后的文件到目标目录
                    final_target_path = os.path.join(dest_dir, target_filename)
                    shutil.move(temp_target_path, final_target_path)
                    print(f"解压至: {final_target_path} 并删除原文件: {gz_file_path}")
                except subprocess.CalledProcessError as e:
                    print(f"解压 {gz_file_path} 时出错: {e}")

file = '/mnt/00.zhibo/random/random1111/'
# 使用示例
decompress_and_move_gz_files(file,file)


