# -*- coding: utf-8 -*-
import warnings

warnings.filterwarnings("ignore")

import os
import glob
import csv
import traceback
import multiprocessing as mp
import MDAnalysis as mda
from MDAnalysis.analysis import align
from spyrmsd import io, rmsd


def find_ref_file(sys_dir):
    """查找参考晶体文件，完美支持 prep.pdb, qp.pdb, *_prep.pdb, *_qp.pdb 等命名"""
    # 1. 优先精确匹配全名文件
    for exact_name in ["prep.pdb", "qp.pdb"]:
        exact_path = os.path.join(sys_dir, exact_name)
        if os.path.exists(exact_path):
            return exact_path

    # 2. 按优先级进行模糊匹配 (依次检查 prep 和 qp)
    prep_patterns = [
        os.path.join(sys_dir, "*_prep.pdb"),
        os.path.join(sys_dir, "*_qp.pdb"),
        os.path.join(sys_dir, "*prep*.pdb"),
        os.path.join(sys_dir, "*qp*.pdb")
    ]
    for pattern in prep_patterns:
        matches = glob.glob(pattern)
        # 过滤掉可能的对接结构文件（如 Structures_*.pdb），防止误吸
        valid_matches = [
            m for m in sorted(matches)
            if not os.path.basename(m).lower().startswith("structure")
        ]
        if valid_matches:
            return valid_matches[0]
    return None


def get_ref_ligand_atoms(ref_universe):
    """提取参考晶体中的配体重原子（排除蛋白、常见溶剂/离子和氢原子）"""
    exclude_resnames = [
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU',
        'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
        'HID', 'HIE', 'HIP', 'CYX', 'ASH', 'GLH', 'LYN', 'MSE', 'ACE', 'NME',
        'HOH', 'WAT', 'DOD', 'SO4', 'PO4', 'EDO', 'GOL', 'DMS', 'ACT', 'CL', 'NA', 'MG', 'ZN', 'CON'
    ]
    ex_str = " ".join(exclude_resnames)
    h_filter = "not name H* and not name [1-9]H* and not element H"

    ligand_atoms = ref_universe.select_atoms(
        f"not protein and not resname {ex_str} and ({h_filter})"
    )
    return ligand_atoms


def get_dock_ligand_atoms(dock_universe, dock_resname="AUTO"):
    """
    智能提取对接文件中的配体重原子
    兼容 MOE ("*")、Covdock ("UNK")、Glide ("LIG")、AutoDock ("UNL") 等各种软件输出
    """
    exclude_resnames = [
        'ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLU', 'GLN', 'GLY', 'HIS', 'ILE', 'LEU',
        'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
        'HID', 'HIE', 'HIP', 'CYX', 'ASH', 'GLH', 'LYN', 'MSE', 'ACE', 'NME',
        'HOH', 'WAT', 'DOD', 'SO4', 'PO4', 'EDO', 'GOL', 'DMS', 'ACT', 'CL', 'NA', 'MG', 'ZN', 'CON'
    ]
    ex_str = " ".join(exclude_resnames)
    h_filter = "not name H* and not name [1-9]H* and not element H"

    # 如果显式指定了有效的残基名且不是通配符
    if dock_resname and dock_resname not in ["AUTO", "*", ""]:
        lig_atoms = dock_universe.select_atoms(f"resname {dock_resname} and ({h_filter})")
        if len(lig_atoms) > 0:
            return lig_atoms

    # 自动识别模式：通过排除蛋白和离子/水分子提取配体
    lig_atoms = dock_universe.select_atoms(
        f"not protein and not resname {ex_str} and ({h_filter})"
    )
    return lig_atoms


def get_matching_protein_heavy_atoms(ref_u, dock_u):
    """
    匹配全部蛋白质重原子（包含 HID, HIE, HIP, CYX 等特有残基，彻底排除氢原子）
    包含 chainID 和 icode 区分，防止多链 (Chain A, B) 或交替构象导致字典 Key 被覆盖
    """
    prot_sel = "protein or resname HID HIE HIP CYX ASH GLH LYN MSE ACE NME"
    h_filter = "not name H* and not name [1-9]H* and not element H"

    full_sel = f"({prot_sel}) and ({h_filter})"

    ref_prot = ref_u.select_atoms(full_sel)
    dock_prot = dock_u.select_atoms(full_sel)

    # 唯一 Key: (链名, 残基号, 插入码, 原子名)
    ref_dict = {(getattr(a, 'chainID', ''), a.resid, getattr(a, 'icode', ''), a.name): a for a in ref_prot}
    dock_dict = {(getattr(a, 'chainID', ''), a.resid, getattr(a, 'icode', ''), a.name): a for a in dock_prot}

    common_keys = sorted(list(set(ref_dict.keys()) & set(dock_dict.keys())))

    if not common_keys:
        return None, None

    ref_align_atoms = mda.AtomGroup([ref_dict[k] for k in common_keys])
    dock_align_atoms = mda.AtomGroup([dock_dict[k] for k in common_keys])

    return ref_align_atoms, dock_align_atoms


def process_single_system(sys_dir, dock_pattern="Structures_*.pdb",
                          dock_resname="AUTO", temp_base_dir="temp_method2_out"):
    sys_name = os.path.basename(os.path.normpath(sys_dir))
    ref_path = find_ref_file(sys_dir)

    # 匹配 Pose 文件并排除 prep/qp/rec/complex 等参考文件
    raw_dock_paths = sorted(glob.glob(os.path.join(sys_dir, dock_pattern)))
    dock_paths = [
        p for p in raw_dock_paths
        if not (p.endswith("prep.pdb") or p.endswith("qp.pdb") or p.endswith("rec.pdb") or p.endswith("complex.pdb"))
    ]
    if not dock_paths:
        dock_paths = sorted(glob.glob(os.path.join(sys_dir, "[Ss]tructure*.pdb")))
        dock_paths = [
            p for p in dock_paths
            if
            not (p.endswith("prep.pdb") or p.endswith("qp.pdb") or p.endswith("rec.pdb") or p.endswith("complex.pdb"))
        ]

    if not ref_path or not os.path.exists(ref_path):
        print(f"[Skip] 未找到晶体参考文件 (prep.pdb / qp.pdb): {sys_dir}")
        return []
    if not dock_paths:
        print(f"[Skip] 未找到对接 Pose 文件 (Structures_*.pdb): {sys_dir}")
        return []

    pid = os.getpid()
    temp_dir = os.path.join(temp_base_dir, f"tmp_{pid}")
    os.makedirs(temp_dir, exist_ok=True)

    system_results = []

    try:
        ref_u = mda.Universe(ref_path)
        ref_lig_atoms = get_ref_ligand_atoms(ref_u)
        n_ref_lig = len(ref_lig_atoms)

        if n_ref_lig == 0:
            print(f"[Warning] 自动提取晶体配体失败: {sys_name} (参考文件: {os.path.basename(ref_path)})")
            return []

        tmp_ref_lig = os.path.join(temp_dir, f"{sys_name}_ref_lig.pdb")
        ref_lig_atoms.write(tmp_ref_lig)
        ref_mol = io.loadmol(tmp_ref_lig)

        for dock_path in dock_paths:
            pose_name = os.path.basename(dock_path).replace(".pdb", "")
            tmp_dock_lig = os.path.join(temp_dir, f"{sys_name}_{pose_name}_aligned_lig.pdb")

            try:
                dock_u = mda.Universe(dock_path)

                # 1. 提取所有重合的蛋白质重原子（除去氢）进行刚体对齐
                ref_align_atoms, dock_align_atoms = get_matching_protein_heavy_atoms(ref_u, dock_u)

                if ref_align_atoms is None or len(ref_align_atoms) == 0:
                    print(f"[Error] {sys_name} -> {pose_name}: 没有匹配的蛋白质重原子")
                    continue

                n_align_atoms = len(ref_align_atoms)

                # 执行刚体对齐
                align.alignto(dock_align_atoms, ref_align_atoms)

                # 2. 提取对接配体重原子（自动识别，兼容 MOE "*" 和 Covdock "UNK"）
                dock_lig_atoms = get_dock_ligand_atoms(dock_u, dock_resname=dock_resname)
                n_dock_lig = len(dock_lig_atoms)

                if n_dock_lig == 0:
                    print(f"[Warning] {sys_name} -> {pose_name}: 未找到配体重原子")
                    continue

                dock_lig_atoms.write(tmp_dock_lig)

                # 3. 计算对齐后的原位对称校正 RMSD
                dock_mol = io.loadmol(tmp_dock_lig)
                val = rmsd.symmrmsd(
                    dock_mol.coordinates,
                    ref_mol.coordinates,
                    dock_mol.atomicnums,
                    ref_mol.atomicnums,
                    dock_mol.adjacency_matrix,
                    ref_mol.adjacency_matrix,
                    minimize=False
                )

                rmsd_val = float(val[0]) if hasattr(val, '__len__') else float(val)

                print(
                    f"[Success] {sys_name} | Ref: {os.path.basename(ref_path)} | Pose: {pose_name} | 对齐蛋白重原子: {n_align_atoms} | 晶体配体重原子: {n_ref_lig} | 对接配体重原子: {n_dock_lig} | RMSD: {round(rmsd_val, 3)} Å")

                system_results.append([
                    sys_name,
                    os.path.basename(ref_path),
                    pose_name,
                    n_align_atoms,
                    n_ref_lig,
                    n_dock_lig,
                    round(rmsd_val, 3)
                ])

            except Exception as e:
                err_detail = str(e) if str(e) else repr(e)
                print(f"[Error] 计算失败 {sys_name} -> {pose_name}: {err_detail}")
            finally:
                if os.path.exists(tmp_dock_lig):
                    os.remove(tmp_dock_lig)

    finally:
        if os.path.exists(tmp_ref_lig):
            os.remove(tmp_ref_lig)

    return system_results


def worker_task(sys_dir):
    return process_single_system(
        sys_dir=sys_dir,
        dock_pattern="Structures_*.pdb",
        dock_resname="AUTO"
    )


if __name__ == "__main__":
    DATASET_DIR = "."
    OUTPUT_CSV = "method2_protein_heavy_aligned_rmsd.csv"

    if find_ref_file(DATASET_DIR):
        system_dirs = [os.path.abspath(DATASET_DIR)]
    else:
        system_dirs = [
            d for d in glob.glob(os.path.join(DATASET_DIR, "*"))
            if os.path.isdir(d) and os.path.basename(d) != "results" and not os.path.basename(d).startswith("temp_")
        ]

    print(f"找到 {len(system_dirs)} 个体系，启动全蛋白重原子对齐 + 对称校正配体 RMSD 计算...")

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
            "Aligned_Protein_Heavy_Atoms",
            "Ref_Ligand_Heavy_Atoms",
            "Dock_Ligand_Heavy_Atoms",
            "Ligand_Symmetry_RMSD_Angstrom"
        ])
        writer.writerows(all_results)

    print(f"\n计算全部完成！共输出 {len(all_results)} 条结果，已保存至 {OUTPUT_CSV}")