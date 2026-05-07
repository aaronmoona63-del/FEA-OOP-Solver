import os
import sys

inp_filename = "/work/tests/Hole_Medium.inp"
out_filename = "/work/tests/mesh_daikonglashen_medium.txt"

if not os.path.exists(inp_filename):
    print(f"❌ 找不到文件 {inp_filename}，请确认已上传！")
    sys.exit(1)

nodes = []
elements = []
reading_nodes = False
reading_elements = False
element_buffer = []

with open(inp_filename, 'r', encoding='utf-8', errors='ignore') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('**'): continue
        upper_line = line.upper()

        if upper_line.startswith('*NODE'):
            reading_nodes = True
            reading_elements = False
            continue
        elif upper_line.startswith('*ELEMENT'):
            if 'C3D8' in upper_line:
                reading_nodes = False
                reading_elements = True
                element_buffer = []
                continue
            else:
                reading_nodes = False
                reading_elements = False
                continue
        elif upper_line.startswith('*'):
            reading_nodes = False
            reading_elements = False
            continue

        if reading_nodes:
            parts = [p.strip() for p in line.split(',')]
            if len(parts) >= 4: nodes.append(parts[:4])

        if reading_elements:
            parts = [p.strip() for p in line.split(',') if p.strip()]
            element_buffer.extend(parts)
            while len(element_buffer) >= 9:
                elements.append(element_buffer[:9])
                element_buffer = element_buffer[9:]

print(f"✅ 解析成功：{len(nodes)} 个节点, {len(elements)} 个单元。")
print(f"💡 适中规模系统总自由度 (DOF): {len(nodes) * 3}")

with open(out_filename, 'w') as f:
    f.write(f"{len(nodes)} {len(elements)}\n")
    for n in nodes: f.write(f"{n[0]} {n[1]} {n[2]} {n[3]}\n")
    for e in elements: f.write(f"{e[0]} {e[1]} {e[2]} {e[3]} {e[4]} {e[5]} {e[6]} {e[7]} {e[8]}\n")
    f.write("210000.0 0.3\n")
