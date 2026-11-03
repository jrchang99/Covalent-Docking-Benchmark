import os
import re


def get_best_energy_run(dlg_file_path):
    runs_data = []
    current_run = None

    # 兼容标准 AD4 与共价对接 (Covalent Docking) 的匹配模式
    run_pattern = re.compile(r'Run\s+(\d+)', re.IGNORECASE)
    # 匹配各类表达：Estimated Free Energy of Binding 或 Binding Energy = -X.XX
    energy_pattern = re.compile(r'(?:Free\s+Energy\s+of\s+Binding|Binding\s+Energy)\s*=\s*([\-\d\.]+)', re.IGNORECASE)

    with open(dlg_file_path, 'r', encoding='utf-8', errors='ignore') as f:
        for line in f:
            # 捕获 Run 编号
            run_match = run_pattern.search(line)
            if run_match and ("Sub-Rank" in line or "USER" in line or "Run" in line):
                current_run = int(run_match.group(1))

            # 捕获自由能
            energy_match = energy_pattern.search(line)
            if energy_match:
                try:
                    energy = float(energy_match.group(1))
                    if current_run is not None:
                        runs_data.append((current_run, energy))
                        current_run = None
                except ValueError:
                    continue

    if not runs_data:
        return None

    # 按自由能升序排序（数值越负，结合力越强，为最佳构象）
    runs_data.sort(key=lambda x: x[1])
    best_run, best_energy = runs_data[0]
    return best_run, best_energy, len(runs_data)


# 递归遍历所有子目录提取
print(f"{'PDB 文件路径':<50} | {'最佳构象(Best Run)':<15} | {'最低自由能 (kcal/mol)':<15}")
print("-" * 85)

for root, dirs, files in os.walk('.'):
    for file in files:
        if file.endswith('.dlg'):
            dlg_path = os.path.join(root, file)
            result = get_best_energy_run(dlg_path)

            rel_path = os.path.relpath(dlg_path, '.')
            if result:
                best_run, best_energy, total_runs = result
                print(f"{rel_path:<50} | Run {best_run:<10} | {best_energy:<15.2f}")
            else:
                print(f"{rel_path:<50} | 未能解析到能量数据")