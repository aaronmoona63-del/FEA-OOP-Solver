cd /work/build

cat << 'EOF' > ../tests/inp2mesh_dense.py
import os

# ==========================================
# ⚙️ 高密网格专属转换脚本 (inp2mesh_dense.py)
# ==========================================

inp_filename = "/work/tests/jiamadaikonglashen.inp"
out_filename = "/work/tests/mesh_daikonglashen_dense.txt"

print(f"正在读取超密 INP 文件: {inp_filename} ...")

nodes = []
elements = []

try:
    with open(inp_filename, 'r') as f:
        lines = f.readlines()
except FileNotFoundError:
    print(f"❌ 找不到文件 {inp_filename}，请确认你已经把它上传到了 /work/tests/ 目录下！")
    exit()

reading_nodes = False
reading_elements = False
element_buffer = []

for line in lines:
    line = line.strip()
    if not line or line.startswith('**'):
        continue
    
    # 判断区块
    if line.startswith('*Node') or line.startswith('*NODE'):
        reading_nodes = True
        reading_elements = False
        continue
    elif line.startswith('*Element') or line.startswith('*ELEMENT'):
        if 'C3D8' in line:  # 确保只读取六面体单元
            reading_nodes = False
            reading_elements = True
            element_buffer = []
            continue
        else:
            reading_nodes = False
            reading_elements = False
            continue
    elif line.startswith('*'):
        # 遇到其他星号开头的配置，停止读取
        reading_nodes = False
        reading_elements = False
        continue
    
    # 提取节点坐标
    if reading_nodes:
        parts = line.split(',')
        if len(parts) >= 4:
            nodes.append([parts[0].strip(), parts[1].strip(), parts[2].strip(), parts[3].strip()])
    
    # 提取单元拓扑 (自动处理 Abaqus 的换行截断)
    if reading_elements:
        parts = [p.strip() for p in line.split(',') if p.strip()]
        element_buffer.extend(parts)
        
        # 1 个单元 ID + 8 个节点 = 9 个数字
        if len(element_buffer) >= 9:
            elements.append(element_buffer[:9])
            element_buffer = []

print(f"✅ 成功解析 {len(nodes)} 个节点, {len(elements)} 个单元。")
print(f"💡 预计系统总自由度 (DOF): {len(nodes) * 3}")

# 写入 OOFEM 标准格式
print(f"正在生成 OOFEM 格式网格文件...")
with open(out_filename, 'w') as f:
    f.write(f"{len(nodes)} {len(elements)}\n")
    
    for n in nodes:
        f.write(f"{n[0]} {n[1]} {n[2]} {n[3]}\n")
        
    for e in elements:
        f.write(f"{e[0]} {e[1]} {e[2]} {e[3]} {e[4]} {e[5]} {e[6]} {e[7]} {e[8]}\n")
        
    # 🌟 终极防报错：在末尾喂给 Fortran 材料属性！
    f.write("210000.0 0.3\n")

print(f"🎉 转换完成！粮草已备齐，文件保存在: {out_filename}")
EOF