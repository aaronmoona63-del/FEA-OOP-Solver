import os
import re
import shutil
import sys

def clean_software_distribution():
    # 目标目录
    target_dir = r"C:\Windows\SoftwareDistribution\Download"
    
    # 匹配正好 32 位的十六进制哈希字符串 (不区分大小写)
    hash_pattern = re.compile(r"^[0-9a-fA-F]{32}$")
    
    if not os.path.exists(target_dir):
        print(f"[错误] 找不到路径: {target_dir}")
        return

    print(f"正在扫描目录: {target_dir}")
    print("-" * 50)

    try:
        # 获取目录下所有条目
        items = os.listdir(target_dir)
        deleted_count = 0

        for item in items:
            item_path = os.path.join(target_dir, item)
            
            # 检查是否为文件夹且符合 32 位哈希正则
            if os.path.isdir(item_path) and hash_pattern.match(item):
                try:
                    print(f"正在删除: {item}")
                    shutil.rmtree(item_path)
                    deleted_count += 1
                except Exception as e:
                    print(f"[失败] 无法删除 {item}: {e}")
            else:
                # 即使是 SharedFileCache，因为不符合正则，也会被安全跳过
                pass

        print("-" * 50)
        print(f"清理完成！共删除了 {deleted_count} 个哈希文件夹。")

    except PermissionError:
        print("[权限拒绝] 请确保以管理员身份运行此脚本！")
    except Exception as e:
        print(f"[未知错误] {e}")

if __name__ == "__main__":
    clean_software_distribution()