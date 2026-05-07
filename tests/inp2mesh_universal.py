import sys, os

if len(sys.argv) < 3:
    print("用法: python inp2mesh_universal.py <输入.inp> <输出.txt>")
    sys.exit(1)

inp_filename = sys.argv[1]
out_filename = sys.argv[2]

nodes, elements = [], []
reading_nodes, reading_elements = False, False
element_buffer = []

with open(inp_filename, 'r', encoding='utf-8', errors='ignore') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('**'): continue
        upper_line = line.upper()

        if upper_line.startswith('*NODE'):
            reading_nodes, reading_elements = True, False
            continue
        elif upper_line.startswith('*ELEMENT'):
            if 'C3D8' in upper_line:
                reading_nodes, reading_elements, element_buffer = False, True, []
                continue
            else:
                reading_nodes, reading_elements = False, False
                continue
        elif upper_line.startswith('*'):
            reading_nodes, reading_elements = False, False
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

with open(out_filename, 'w') as f:
    f.write(f"{len(nodes)} {len(elements)}\n")
    for n in nodes: f.write(f"{n[0]} {n[1]} {n[2]} {n[3]}\n")
    for e in elements: f.write(f"{e[0]} {e[1]} {e[2]} {e[3]} {e[4]} {e[5]} {e[6]} {e[7]} {e[8]}\n")
    f.write("210000.0 0.3\n")

print(f"✅ {inp_filename} 解析成功: {len(nodes)} 节点, 自由度 {len(nodes)*3}")
