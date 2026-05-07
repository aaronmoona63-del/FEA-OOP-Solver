# =======================================================
# ⚙️ 60万自由度均匀板网格生成器 (严格修复版)
# =======================================================

nx, ny, nz = 210, 105, 8
lx, ly, lz = 200.0, 100.0, 5.0
out_filename = "/work/tests/mesh_regular_dense.txt"

print(f">>> 正在生成超大规模均匀网格 (严格匹配 Fortran 读取格式)...")
nodes = []
for k in range(nz + 1):
    z = k * (lz / nz)
    for j in range(ny + 1):
        y = -50.0 + j * (ly / ny)
        for i in range(nx + 1):
            x = -100.0 + i * (lx / nx)
            nodes.append((x, y, z))

def get_node_id(i, j, k):
    return k * (nx + 1) * (ny + 1) + j * (nx + 1) + i + 1

elements = []
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
            elements.append((n1, n2, n3, n4, n5, n6, n7, n8))

nnodes = len(nodes)
nelems = len(elements)

with open(out_filename, 'w') as f:
    f.write(f"{nnodes} {nelems}\n")
    # 🌟 核心修复 1：加入 idx+1 作为节点 ID，补齐 4 列
    for idx, n in enumerate(nodes):
        f.write(f"{idx+1} {n[0]:.6f} {n[1]:.6f} {n[2]:.6f}\n")
        
    # 🌟 核心修复 2：加入 idx+1 作为单元 ID，补齐 9 列
    for idx, e in enumerate(elements):
        f.write(f"{idx+1} {e[0]} {e[1]} {e[2]} {e[3]} {e[4]} {e[5]} {e[6]} {e[7]}\n")
        
    # 补充材料属性
    f.write("210000.0 0.3\n")

print(f"✅ 完美格式网格已生成！节点数: {nnodes}, 单元数: {nelems}")
