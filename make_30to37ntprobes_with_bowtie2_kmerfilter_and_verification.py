# -*- coding: utf-8 -*-
import subprocess
import sys
import os
import pandas as pd
import numpy as np
import shutil
import argparse # argparse をインポート

# === 追加: 逆相補鎖を生成するヘルパー関数 ===
def get_reverse_complement(dna_seq):
    """DNA配列の逆相補鎖を返す"""
    complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    try:
        # 大文字に統一し、マップを使って相補鎖を作成
        complement_seq = "".join(complement_map[base.upper()] for base in dna_seq)
        # 文字列を逆順にして返す
        return complement_seq[::-1]
    except KeyError as e:
        print(f"  [検証エラー] 不正な塩基が含まれています: {e}")
        return None
# === 追加ここまで ===


# --- 引数の設定 ---
# sys.argv を直接参照する代わりに argparse を使用
parser = argparse.ArgumentParser(
    description="Probe generation pipeline with bowtie2 alignment and k-mer filtering.",
    usage="%(prog)s [input_folder] [bowtie2_index] [jellyfish_file] [-m MER_LENGTH] [-k KMER_THRESHOLD]"
)

# 必須の「位置引数」
parser.add_argument(
    "input_folder", 
    help="Input folder containing target .txt files."
)
parser.add_argument(
    "bowtie2_index", 
    help="Path to the bowtie2 index prefix (e.g., /path/to/hg38)."
)
parser.add_argument(
    "jellyfish_file", 
    help="Path to the Jellyfish .jf file for kmerFilter (e.g., /path/to/hg38_20mer.jf)."
)

# オプション引数 (kmerFilter用)
parser.add_argument(
    "-m", "--mer_length", 
    default="18", # kmerFilter_py3.py のデフォルト値に合わせる
    help="The k-mer length (merLength) for kmerFilter. (Default: 18)"
)
parser.add_argument(
    "-k", "--kmer_threshold", 
    default="5",  # kmerFilter_py3.py のデフォルト値に合わせる
    help="The k-mer occurrence threshold (kmerThreshold) for kmerFilter. (Default: 5)"
)

# 引数のパース実行
if len(sys.argv) == 1:
    parser.print_help(sys.stderr)
    sys.exit(1)

args = parser.parse_args()

# パースした引数を変数に割り当て
input_file_folder = args.input_folder
bowtie2_index_path = args.bowtie2_index
jellyfish_file_path = args.jellyfish_file
mer_length = args.mer_length
kmer_threshold = args.kmer_threshold

print(f"--- パイプライン開始 ---")
print(f"  Input folder: {input_file_folder}")
print(f"  Bowtie2 index: {bowtie2_index_path}")
print(f"  Jellyfish file: {jellyfish_file_path}")
print(f"  k-mer length (m): {mer_length}")
print(f"  k-mer threshold (k): {kmer_threshold}")
print("------------------------")

# --- 出力フォルダの設定 ---
fastq_folder_path = os.path.join(input_file_folder, 'fastq')
sam_folder_path = os.path.join(input_file_folder, 'sam')
bed_folder_path = os.path.join(input_file_folder, 'bed')
kfiltered_bed_folder_path = os.path.join(input_file_folder, 'kfiltered_bed')
upto20_folder_path = os.path.join(input_file_folder, 'upto20')

# --- 出力フォルダの作成 ---
if not os.path.exists(fastq_folder_path):
    os.makedirs(fastq_folder_path)
if not os.path.exists(sam_folder_path):
    os.makedirs(sam_folder_path)
if not os.path.exists(bed_folder_path):
    os.makedirs(bed_folder_path)
if not os.path.exists(kfiltered_bed_folder_path):
    os.makedirs(kfiltered_bed_folder_path)
if not os.path.exists(upto20_folder_path):
    os.makedirs(upto20_folder_path)

# "targets" フォルダ内の .txt ファイルを取得
txt_files = [f for f in os.listdir(input_file_folder) if os.path.isfile(os.path.join(input_file_folder, f)) and f.endswith('.txt')]

# --- メイン処理ループ ---
for txt_file in txt_files:
    print(f"\n--- 処理中: {txt_file} ---")
    txt_path = os.path.join(input_file_folder, txt_file)
    
    # 1. blockParse.py の実行
    print(f"  [1/6] blockParse を実行中...") # ステップ数を 1/6 に変更
    fastq_file = txt_file.replace('.txt', '')
    fastq_path = os.path.join(fastq_folder_path, fastq_file)
    subprocess.call(['python', 'blockParse_py3.py','-l', '30','-L', '37','-f', txt_path,'-o', fastq_path])
    
    # 2. bowtie2 の実行
    print(f"  [2/6] bowtie2 を実行中...") # ステップ数を 2/6 に変更
    fastq_file_ext = txt_file.replace('.txt', '.fastq')
    fastq_path_ext = os.path.join(fastq_folder_path, fastq_file_ext)
    sam_file = txt_file.replace('.txt', '.sam')
    sam_path = os.path.join(sam_folder_path, sam_file)
    subprocess.call(['bowtie2', '-x', bowtie2_index_path, '-U', fastq_path_ext, '--no-hd','-t', '-k', '100', '-D', '25', '-R', '3', '-N', '0', '-L', '20', '-i', 'S,1,0.50', '-S', sam_path])
    
    # 3. OutputClean.py の実行
    print(f"  [3/6] OutputClean を実行中...") # ステップ数を 3/6 に変更
    bed_file = txt_file.replace('.txt', '')
    bed_path = os.path.join(bed_folder_path, bed_file)
    subprocess.call(['python', 'outputClean_py3_ver2.py', '-uz','-f', sam_path, '-o', bed_path])

    # OutputClean が出力した .bed ファイルのフルパス
    bed_file_ext = txt_file.replace('.txt', '.bed')
    bed_path_ext = os.path.join(bed_folder_path, bed_file_ext)

    # 4. kmerFilter_py3.py の実行
    print(f"  [4/6] kmerFilter を実行中 (k={kmer_threshold}, m={mer_length})...") # ステップ数を 4/6 に変更
    kfiltered_bed_file = txt_file.replace('.txt', '_kfiltered')
    kfiltered_bed_path_prefix = os.path.join(kfiltered_bed_folder_path, kfiltered_bed_file)
    
    subprocess.call([
        'python', 'kmerFilter_py3.py',
        '-f', bed_path_ext,
        '-j', jellyfish_file_path,
        '-m', mer_length,
        '-k', kmer_threshold,
        '-o', kfiltered_bed_path_prefix
    ])

    # kmerFilter が出力した .bed ファイルのフルパス
    kfiltered_bed_path_ext = f"{kfiltered_bed_path_prefix}.bed"

    # 5. upto20 (プローブ数の絞り込み) の実行
    print(f"  [5/6] upto20 (プローブ数絞り込み) を実行中...") # ステップ数を 5/6 に変更
    
    try:
        df = pd.read_csv(kfiltered_bed_path_ext, sep='\t', header=None)
    except pd.errors.EmptyDataError:
        print(f"  警告: kmerFilter の結果、{kfiltered_bed_path_ext} が空になりました。upto20処理をスキップします。")
        continue 
    except FileNotFoundError:
        print(f"  エラー: kmerFilter が {kfiltered_bed_path_ext} を生成しませんでした。OutputCleanの結果が空だった可能性があります。スキップします。")
        continue

    upto20_file = txt_file.replace('.txt', '_RCupto20.bed')
    upto20_path = os.path.join(upto20_folder_path, upto20_file)

    if len(df) <= 20:
        shutil.copyfile(kfiltered_bed_path_ext, upto20_path)
    else:
        df['length'] = df[3].apply(lambda x: len(str(x)))
        median = np.median(df['length'])
        df['diff'] = abs(df['length'] - median)
        df = df.nsmallest(20, 'diff')
        df = df.sort_values(by=[1])
        df.to_csv(upto20_path, sep='\t', header=False, index=False, columns=[0, 1, 2, 3, 4])

    
    # === [6/6] 検証ステップを追加 ===
    print(f"  [6/6] 最終プローブの検証を実行中...")
    
    # 1. 元の配列を .txt から読み込む
    target_seq = ""
    try:
        with open(txt_path, 'r') as f_target:
            for line in f_target:
                # ヘッダー行 ('>') をスキップ
                if not line.startswith('>'):
                    target_seq += line.strip().upper()
    except FileNotFoundError:
        print(f"  [検証エラー] 元の .txt ファイルが見つかりません: {txt_path}")
        continue # このファイルの処理を中断し、次の .txt ファイルへ

    # 2. 生成されたプローブを _RCupto20.bed から読み込む
    probe_sequences = []
    try:
        with open(upto20_path, 'r') as f_probes:
            for line in f_probes:
                parts = line.strip().split('\t')
                if len(parts) >= 4:
                    probe_sequences.append(parts[3].upper())
    except FileNotFoundError:
        # kmerFilter等で0個になった場合はファイルが生成されないので、これはエラーではない
        print(f"  [検証スキップ] プローブファイルが生成されていません: {upto20_path}")
        print(f"--- 完了 (プローブ0個): {txt_file} ---")
        continue # 次のファイルの処理へ

    # 3. 検証を実行
    total_probes = len(probe_sequences)
    rc_found_count = 0
    
    if total_probes > 0 and target_seq:
        for probe_seq in probe_sequences:
            # プローブの「逆相補鎖」を計算
            probe_rc = get_reverse_complement(probe_seq)
            
            # 逆相補鎖が「元の配列」に含まれているか確認
            if probe_rc and probe_rc in target_seq:
                rc_found_count += 1
            else:
                # 含まれていない場合（デバッグ用にNGのプローブを表示）
                print(f"    [検証NG] プローブ: {probe_seq} (RC: {probe_rc}) が元の配列に見つかりません。")

        # 4. 検証結果をレポート
        print(f"  [検証結果] {rc_found_count} / {total_probes} 個のプローブが、元の配列の逆相補鎖として見つかりました。")
        if rc_found_count != total_probes:
            print(f"  [警告] {total_probes - rc_found_count} 個のプローブが逆相補鎖として見つかりませんでした。")
            
    elif total_probes == 0:
        print("  [検証スキップ] 生成されたプローブが0個でした。")
    else:
        # target_seq が空だった場合
        print("  [検証エラー] 元の配列の読み込みに失敗しました（空でした）。")
    # === 検証ロジックここまで ===

    print(f"--- 完了: {txt_file} ---")

print("\nすべての処理が完了しました。")