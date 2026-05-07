import numpy as np

# ==========================================
# ⚙️ 规则矩形板网格生成器 (修复 EOF 报错版)
# ==========================================

# 1. 物理尺寸定义
Length = 200.0   
Width = 100.0    
Thickness = 5.0  

# 2. 网格划分密度
nx = 50 
ny = 40 
nz = 4   

# 3. 计算节点与单元总数
nnodes_x, nnodes_y, nnodes_z = nx + 1, ny + 1, nz + 1
total_nodes = nnodes_x * nnodes_y * nnodes_z
total_elems = nx * ny * nz

print(f"正在生成完美规则网格...")
print(f"单元数: {total_elems}")
print(f"节点数: {total_nodes} (自由度: {total_nodes * 3})")

# 4. 生成节点坐标
x_coords = np.linspace(-100.0, 100.0, nnodes_x)
y_coords = np.linspace(-50.0, 50.0, nnodes_y)
z_coords = np.linspace(0.0, 5.0, nnodes_z)

def get_node_id(i, j, k):
    return k * (nnodes_x * nnodes_y) + j * nnodes_x + i + 1

# 5. 写入 TXT 文件
filename = "mesh_regular_plate.txt"
with open(filename, "w") as f:
    # 写入文件头
    f.write(f"{total_nodes} {total_elems}\n")
    
    # 写入所有节点
    for k in range(nnodes_z):
        for j in range(nnodes_y):
            for i in range(nnodes_x):
                node_id = get_node_id(i, j, k)
                f.write(f"{node_id} {x_coords[i]:.6f} {y_coords[j]:.6f} {z_coords[k]:.6f}\n")
                
    # 写入所有单元
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
                
    # 🌟 核心修复：在文件最末尾写入材料属性 (E 和 nu)，满足 Fortran 的读取胃口！
    f.write("210000.0 0.3\n")

print(f"✅ 生成完毕！文件已保存为: {filename}")