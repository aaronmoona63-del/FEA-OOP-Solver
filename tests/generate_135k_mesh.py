import sys

def generate_mesh(nx, ny, nz, lx, ly, lz, filename, E, nu, include_ids=True):
    n_nodes = (nx + 1) * (ny + 1) * (nz + 1)
    n_elems = nx * ny * nz
    print(f"正在生成网格：{n_elems} 个单元, {n_nodes} 个节点...")
    
    nodes_x = nx + 1
    nodes_y = ny + 1
    nodes_z = nz + 1
    
    try:
        with open(filename, 'w') as f:
            # 1. 写入 Header: 节点数 单元数
            f.write(f"{n_nodes} {n_elems}\n")
            
            # 2. 写入所有节点坐标 (带有节点编号)
            print("正在写入节点坐标...")
            node_id = 1
            for k in range(nodes_z):
                z = k * lz / nz
                for j in range(nodes_y):
                    y = j * ly / ny
                    for i in range(nodes_x):
                        x = i * lx / nx
                        if include_ids:
                            # 格式：编号 X Y Z
                            f.write(f"{node_id} {x:.6f} {y:.6f} {z:.6f}\n")
                        else:
                            f.write(f"{x:.6f} {y:.6f} {z:.6f}\n")
                        node_id += 1
                        
            # 3. 写入所有单元连接关系 (带有单元编号)
            print("正在写入单元拓扑关系...")
            elem_id = 1
            for k in range(nz):
                for j in range(ny):
                    for i in range(nx):
                        n1 = i + j * nodes_x + k * nodes_x * nodes_y + 1
                        n2 = n1 + 1
                        n3 = n2 + nodes_x
                        n4 = n1 + nodes_x
                        n5 = n1 + nodes_x * nodes_y
                        n6 = n2 + nodes_x * nodes_y
                        n7 = n3 + nodes_x * nodes_y
                        n8 = n4 + nodes_x * nodes_y
                        if include_ids:
                            # 格式：编号 N1 N2 ... N8
                            f.write(f"{elem_id} {n1} {n2} {n3} {n4} {n5} {n6} {n7} {n8}\n")
                        else:
                            f.write(f"{n1} {n2} {n3} {n4} {n5} {n6} {n7} {n8}\n")
                        elem_id += 1
            
            # 4. 追加材料属性 (杨氏模量 泊松比)
            print(f"正在追加材料属性: E = {E}, nu = {nu} ...")
            f.write(f"{E:.6f} {nu:.6f}\n")
                        
        print(f"✅ 格式修正完毕！文件已成功保存为：{filename}")
    except Exception as e:
        print(f"❌ 发生错误：{e}")

if __name__ == '__main__':
    # 调用函数，注意 include_ids=True 开启了编号写入！
    generate_mesh(nx=30, ny=30, nz=150, 
                  lx=10.0, ly=10.0, lz=50.0, 
                  filename='/work/tests/mesh_scale.txt',
                  E=210000.0, nu=0.3,
                  include_ids=True)