import os
import glob

def convert_inp_to_txt():
    print("="*55)
    print(" 🚀 Abaqus INP -> Fortran OOP 网格转换器 🚀")
    print("="*55)

    # 🌟 绝杀功能 1：自动寻找当前目录下的所有 .inp 文件
    inp_files = glob.glob("*.inp")
    
    if not inp_files:
        print("❌ 致命错误：在当前目录下没有找到任何 .inp 文件！")
        print("💡 抢救指南：请确认你已经把 Abaqus 导出的文件放在了与本脚本相同的文件夹下。")
        return

    # 默认处理找到的第一个 .inp 文件
    inp_file = inp_files[0]
    if len(inp_files) > 1:
        print(f"⚠️ 发现多个 .inp 文件，将默认处理第一个: {inp_file}")
    else:
        print(f"🔍 成功锁定 Abaqus 网格文件: {inp_file}")

    # 🌟 绝杀功能 2：自动生成配套的输出文件名 (比如 Job-HolePlate.inp -> mesh_Job-HolePlate.txt)
    base_name = os.path.splitext(inp_file)[0]
    out_file = f"mesh_{base_name}.txt"

    nodes = []
    elements = []

    print("⏳ 正在高速解析底层拓扑数据，请稍候...")

    # 兼容处理各种编码格式
    try:
        with open(inp_file, 'r', encoding='utf-8') as f:
            lines = f.readlines()
    except UnicodeDecodeError:
        with open(inp_file, 'r', encoding='gbk') as f:
            lines = f.readlines()

    state = 0 # 状态机：0-搜索中, 1-读取节点, 2-读取单元

    for line in lines:
        line = line.strip()
        
        # 跳过空行和 Abaqus 的注释行
        if not line or line.startswith('**'):
            continue
        
        # 捕捉节点块和单元块 (不区分大小写，更鲁棒)
        if line.upper().startswith('*NODE'):
            state = 1
            continue
        elif line.upper().startswith('*ELEMENT'):
            if 'C3D8' not in line.upper():
                print("⚠️ 警告：检测到的单元可能不是 C3D8 (Hex8) 类型，可能会导致 Fortran 组装失败！")
            state = 2
            continue
        elif line.startswith('*'): 
            # 遇到其他关键字组合（如 *NSET, *MATERIAL），立刻退出读取状态
            state = 0
            continue

        # 🌟 绝杀功能 3：更鲁棒的字符串解析（防范 Abaqus 偶尔在行尾多加逗号的 Bug）
        if state == 1:
            parts = [p.strip() for p in line.split(',')]
            parts = [p for p in parts if p] # 过滤掉空字符串
            if len(parts) >= 4:
                nodes.append(f"{parts[0]} {parts[1]} {parts[2]} {parts[3]}")
                
        elif state == 2:
            parts = [p.strip() for p in line.split(',')]
            parts = [p for p in parts if p] # 过滤掉空字符串
            # Hex8 单元必须包含：ID, 加上 8个节点编号 (共9个数据)
            if len(parts) >= 9:
                # 只取前9个，防止有些软件导出时附带额外信息
                elements.append(" ".join(parts[:9]))

    if len(nodes) == 0 or len(elements) == 0:
        print("❌ 解析失败：未能从文件中提取到有效的节点或单元数据，请检查 INP 文件格式！")
        return

    # 写入 Fortran 求解器专属的 txt 格式
    print(f"💾 正在将规范化数据封装至 {out_file} ...")
    with open(out_file, 'w', encoding='utf-8') as f:
        # 第一行：节点总数 单元总数
        f.write(f"{len(nodes)} {len(elements)}\n")
        
        # 写入节点坐标
        for n in nodes:
            f.write(n + "\n")
        
        # 写入单元拓扑
        for e in elements:
            f.write(e + "\n")
            
        # 写入默认材料属性 (钢材: 弹性模量 E=210000, 泊松比 nu=0.3)
        f.write("210000.0 0.3\n")

    print("\n" + "="*55)
    print(" 🎉 格式转换完美竣工！")
    print(f" 📊 数据统计:")
    print(f"    - 有效节点数 (Nodes): {len(nodes):,}")
    print(f"    - Hex8单元数 (Elems): {len(elements):,}")
    print(f" 🎯 Fortran 挂载目标文件: {out_file}")
    print("="*55 + "\n")

if __name__ == "__main__":
    convert_inp_to_txt()