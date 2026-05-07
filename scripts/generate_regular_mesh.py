import numpy as np
import os

# ==========================================
# ⚙️ 规则矩形板网格生成器 (强制绝对路径版)
# ==========================================

Length = 200.0   
Width = 100.0    
Thickness = 5.0  

nx = 50 
ny = 40 
nz = 4   

nnodes_x, nnodes_y, nnodes_z = nx + 1, ny + 1, nz + 1
total_nodes = nnodes_x * nnodes_y * nnodes_z
total_elems = nx * ny * nz

print(f"正在生成完美规则网格...")
print(f"单元数: {total_elems}")
print(f"节点数: {total_nodes} (自由度: {total_nodes * 3})")

x_coords = np.linspace(-100.0, 100.0, nnodes_x)
y_coords = np.linspace(-50.0, 50.0, nnodes_y)
z_coords = np.linspace(0.0, 5.0, nnodes_z)

def get_node_id(i, j, k):
    return k * (nnodes_x * nnodes_y) + j * nnodes_x + i + 1

# 🌟 核心修复：强制将文件生成到 Fortran 指定读取的 tests 目录下！
filename = "/work/tests/mesh_regular_plate.txt"

# 确保 tests 文件夹存在
os.makedirs(os.path.dirname(filename), exist_ok=True)

with open(filename, "w") as f:
    f.write(f"{total_nodes} {total_elems}\n")
    
    for k in range(nnodes_z):
        for j in range(nnodes_y):
            for i in range(nnodes_x):
                node_id = get_node_id(i, j, k)
                f.write(f"{node_id} {x_coords[i]:.6f} {y_coords[j]:.6f} {z_coords[k]:.6f}\n")
                
    elem_id = 1
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                n1 = get_node_id(i, j, k)
                n2 = get_node_id(i+1, j, k)
                n3 = get_node_id(i+1, j+1, k)
                n4 = get_node_id(i, j+1, k)
                n5 = get_node_id(i, j, k+1)
                n6 = get_node_id(i+1, j, k+1)
                n7 = get_node_id(i+1, j+1, k+1)
                n8 = get_node_id(i, j+1, k+1)
                f.write(f"{elem_id} {n1} {n2} {n3} {n4} {n5} {n6} {n7} {n8}\n")
                elem_id += 1
                
    # 写入材料属性，满足 Fortran 读取胃口
    f.write("210000.0 0.3\n")

print(f"✅ 生成完毕！文件已强行覆盖保存至: {filename}")