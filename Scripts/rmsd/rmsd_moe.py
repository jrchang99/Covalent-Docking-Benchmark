# -*- coding: utf-8 -*-
import warnings

warnings.filterwarnings("ignore")

import os
import glob
import csv
import traceback
import multiprocessing as mp
import numpy as np
import MDAnalysis as mda
from spyrmsd import io, rmsd


def find_ref_file(sys_dir):
    """优先匹配 lig_ref.pdb"""
    exact_path = os.path.join(sys_dir, "lig_ref.pdb")
    if os.path.exists(exact_path):
        return exact_path

    for exact_name in ["ref.pdb", "prep.pdb", "qp.pdb"]:
        p = os.path.join(sys_dir, exact_name)
        if os.path.exists(p):
            return p

    patterns = [
        os.path.join(sys_dir, "*ref*.pdb"),
        os.path.join(sys_dir, "*prep*.pdb"),
        os.path.join(sys_dir, "*qp*.pdb")
    ]
    for pattern in patterns:
        matches = glob.glob(pattern)
        valid = [
            m for m in sorted(matches)
            if not os.path.basename(m).lower().startswith("structure")
        ]
        if valid:
            return valid[0]
    return None


def get_ref_ligand_atoms(ref_u):
    """提取 lig_ref.pdb 中的纯配体重原子 (HETATM / 非蛋白)"""
    h_filter = "not name H* and not name [1-9]H* and not element H"
    ref_lig_atoms = ref_u.select_atoms(f"not protein and record_type HETATM and ({h_filter})")
    if len(ref_lig_atoms) == 0:
        ref_lig_atoms = ref_u.select_atoms(f"not record_type ATOM and ({h_filter})")
    return ref_lig_atoms


def extract_moe_dock_ligand(ref_lig_atoms, dock_u):
    """
    根据 MOE 文件的结构特征精准提取纯配体重原子：
    MOE 对接文件的最后 N_lig 个重原子必定是配体。
    """
    h_filter = "not name H* and not name [1-9]H* and not element H"
    dock_heavy = dock_u.select_atoms(h_filter)

    n_ref_lig = len(ref_lig_atoms)
    n_dock_heavy = len(dock_heavy)

    if n_dock_heavy < n_ref_lig:
        print(f"[Warning] MOE 文件重原子总数 ({n_dock_heavy}) 小于 参考配体重原子数 ({n_ref_lig})")
        return None

    # 取 MOE 文件中最后 n_ref_lig 个重原子（纯配体部分）
    dock_lig_atoms = dock_heavy[-n_ref_lig:]
    return dock_lig_atoms


def load_spyrmsd_mol(path):
    """加载分子并兼容 spyrmsd 返回单个 Molecule 或 List 的情况"""
    mol = io.loadmol(path)
    if isinstance(mol, list):
        return mol[0]
    return mol


def process_single_system(sys_dir):
    sys_name = os.path.basename(os.path.normpath(sys_dir))
    ref_path = find_ref_file(sys_dir)

    dock_paths = sorted(glob.glob(os.path.join(sys_dir, "[Ss]tructure*.pdb")))
    if not dock_paths:
        raw_dock_paths = sorted(glob.glob(os.path.join(sys_dir, "*.pdb")))
        dock_paths = [
            p for p in raw_dock_paths
            if not any(k in os.path.basename(p).lower() for k in ["ref", "prep", "qp", "rec", "complex"])
        ]

    if not ref_path or not os.path.exists(ref_path):
        print(f"[Skip] 未找到参考文件 lig_ref.pdb: {sys_dir}")
        return []
    if not dock_paths:
        print(f"[Skip] 未找到对接 Pose 文件 (structures_*.pdb): {sys_dir}")
        return []

    system_results = []

    try:
        ref_u = mda.Universe(ref_path)
        ref_lig_atoms = get_ref_ligand_atoms(ref_u)

        if len(ref_lig_atoms) == 0:
            print(f"[Error] 参考文件 lig_ref.pdb 未提取到配体重原子 | {sys_name}")
            return []

        # 1. 从原始 lig_ref.pdb 直接加载 spyrmsd 分子结构与拓扑信息
        spy_ref_mol = load_spyrmsd_mol(ref_path)
        ref_lig_indices = [a.index for a in ref_lig_atoms]

        ref_coords = ref_lig_atoms.positions
        ref_atomicnums = spy_ref_mol.atomicnums[ref_lig_indices]
        ref_adj = spy_ref_mol.adjacency_matrix[np.ix_(ref_lig_indices, ref_lig_indices)]

        for dock_path in dock_paths:
            pose_name = os.path.basename(dock_path).replace(".pdb", "")

            try:
                dock_u = mda.Universe(dock_path)

                # 2. 从 MOE 对接文件中提取纯配体重原子
                dock_lig_atoms = extract_moe_dock_ligand(ref_lig_atoms, dock_u)

                if dock_lig_atoms is None or len(dock_lig_atoms) != len(ref_lig_atoms):
                    print(
                        f"[Error] {sys_name} -> {pose_name}: 配体重原子数不一致 (Dock={len(dock_lig_atoms) if dock_lig_atoms else 0}, Ref={len(ref_lig_atoms)})")
                    continue

                # 3. 直接获取 MOE 原始坐标（不对坐标做任何平移或旋转对齐）
                dock_coords = dock_lig_atoms.positions

                # 4. 从原始 structures_*.pdb 加载分子拓扑信息
                spy_dock_mol = load_spyrmsd_mol(dock_path)
                dock_lig_indices = [a.index for a in dock_lig_atoms]

                dock_atomicnums = spy_dock_mol.atomicnums[dock_lig_indices]
                dock_adj = spy_dock_mol.adjacency_matrix[np.ix_(dock_lig_indices, dock_lig_indices)]

                # 5. 计算原位 (In-place) 对称校正配体 RMSD
                val = rmsd.symmrmsd(
                    dock_coords,
                    ref_coords,
                    dock_atomicnums,
                    ref_atomicnums,
                    dock_adj,
                    ref_adj,
                    minimize=False  # minimize=False 表示不对齐坐标，直接计算原始位置 RMSD
                )

                ligand_rmsd_val = float(val[0]) if hasattr(val, '__len__') else float(val)

                print(
                    f"[Success] {sys_name} | Pose: {pose_name} | "
                    f"晶体配体重原子: {len(ref_lig_atoms)} | 对接配体重原子: {len(dock_lig_atoms)} | "
                    f"配体原位对称 RMSD: {round(ligand_rmsd_val, 3)} Å"
                )

                system_results.append([
                    sys_name,
                    os.path.basename(ref_path),
                    pose_name,
                    len(ref_lig_atoms),
                    len(dock_lig_atoms),
                    round(ligand_rmsd_val, 3)
                ])

            except Exception as e:
                err_detail = str(e) if str(e) else repr(e)
                print(f"[Error] 计算失败 {sys_name} -> {pose_name}: {err_detail}")
                traceback.print_exc()

    except Exception as e:
        err_detail = str(e) if str(e) else repr(e)
        print(f"[Error] 处理体系失败 {sys_name}: {err_detail}")
        traceback.print_exc()

    return system_results


def worker_task(sys_dir):
    return process_single_system(sys_dir)


if __name__ == "__main__":
    DATASET_DIR = "."
    OUTPUT_CSV = "moe_ligand_only_inplace_rmsd.csv"

    if find_ref_file(DATASET_DIR):
        system_dirs = [os.path.abspath(DATASET_DIR)]
    else:
        system_dirs = [
            d for d in glob.glob(os.path.join(DATASET_DIR, "*"))
            if os.path.isdir(d) and not os.path.basename(d).startswith("temp_")
        ]

    print(f"找到 {len(system_dirs)} 个体系，开始执行纯配体原位对称 RMSD 计算 (无对齐)...")

    num_processes = min(mp.cpu_count(), 16)
    all_results = []

    with mp.Pool(processes=num_processes) as pool:
        for res in pool.imap_unordered(worker_task, system_dirs):
            if res:
                all_results.extend(res)

    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "System",
            "Ref_File",
            "Pose",
            "Ref_Ligand_Heavy_Atoms",
            "Dock_Ligand_Heavy_Atoms",
            "Ligand_InPlace_Symmetry_RMSD_Angstrom"
        ])
        writer.writerows(all_results)

    print(f"\n全部计算完成！输出结果已保存至 {OUTPUT_CSV}")