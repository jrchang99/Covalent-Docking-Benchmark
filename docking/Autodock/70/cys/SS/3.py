import os
import re


def get_highest_energy_pose(dlg_file):
    run_energies = []
    # 匹配 AutoDock dlg 中的 Binding Energy 行
    energy_pattern = re.compile(r'RANKING\s+(\d+)\s+(\d+)\s+([\d\.\-]+)')

    with open(dlg_file, 'r', errors='ignore') as f:
        for line in f:
            # 方式1：匹配 RANKING 表格 (Rank, Run, Binding Energy)
            match = energy_pattern.search(line)
            if match:
                rank, run, energy = match.groups()
                run_energies.append((int(run), float(energy)))

            # 方式2：备用匹配模式
            elif "Estimated Free Energy of Binding" in line:
                parts = line.split("=")
                if len(parts) > 1:
                    try:
                        e_val = float(parts[1].split()[0])
                        # 提取 Run 编号
                        # 此处配合文件上下文使用
                    except ValueError:
                        pass

    if run_energies:
        # 按能量从大到小排序 (能量数值越大，代表能量越高/结合越差)
        run_energies.sort(key=lambda x: x[1], reverse=True)
        highest_run, highest_e = run_energies[0]
        print(f"文件: {os.path.basename(dlg_file)} --> 最高能量结构为: Run {highest_run} (能量: {highest_e} kcal/mol)")
    else:
        print(f"未在 {os.path.basename(dlg_file)} 中匹配到能量信息")


# 示例：检测当前目录下的 dlg 文件
for root, dirs, files in os.walk("."):
    for file in files:
        if file.endswith(".dlg"):
            get_highest_energy_pose(os.path.join(root, file))