from Bio import PDB

def PDB_CAP(input_dir, output_dir):
    """
    PDBファイルを読み込み、N末端をACE、C末端をNMEに変更する関数
    :param input_dir: 入力PDBファイルのパス
    :param output_dir: 出力PDBファイルのパス
    """

    # PDBファイルの読み込み
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure('Protein', input_dir)

    dresseq = 0.1 #アミノ酸残基番号の変位(デフォルト0変位だと残基番号の被りが発生するため)

    ACEAtom, NMEAtom = None, None #N末、C末の水素原子を格納する変数
    ACEres = None #ACE残基を格納する変数
    NMEres = None #NME残基を格納する変数

    for model in structure:
        for chain in model:
            for residue in chain:
                hasOXT = False #OXT原子の有無を確認するフラグ
                atoms_to_delete = []
                for atom in residue:
                    if atom.get_name().strip() == 'HN1':
                        # N末水素原子の名前を変更
                        ACEAtom = atom # ACE残基に追加する水素原子
                        atoms_to_delete.append(atom.id) # 元の残基から削除
                        dresseq += 1

                        ACEres_id = residue.get_id()
                        ACEres = PDB.Residue.Residue((' ', ACEres_id[1] + dresseq - 1, ' '), 'ACE', ' ')
                        ACEAtom.fullname = 'C'
                        ACEres.add(ACEAtom)
                        
                    if (atom.get_name().strip() == 'HN2') or (atom.get_name().strip() == 'HN3'):
                        atoms_to_delete.append(atom.id) # 元の残基から削除

                    if atom.get_name().strip() == 'OXT':
                        NMEAtom = atom #NME残基に追加する酸素原子
                        atoms_to_delete.append(atom.id)
                        dresseq += 0

                        NMEres_id = residue.get_id()
                        NMEres = PDB.Residue.Residue((' ', NMEres_id[1] + dresseq + 1, ' '), 'NME', ' ')
                        NMEAtom.fullname = 'N'
                        NMEres.add(NMEAtom)
                        hasOXT = True
                
                #不要な原子を一括削除
                for atom_id in atoms_to_delete:
                    residue.detach_child(atom_id)
                
                newresidue_id = residue.get_id()
                newresidue_id = (newresidue_id[0], newresidue_id[1] + dresseq, newresidue_id[2])
                residue.id = newresidue_id #残基番号の変更

                dresseq += hasOXT #OXT原子があった場合は次から残基番号を1増加

            if ACEres is not None:
                chain.add(ACEres)
                ACEres = None
            
            if NMEres is not None:
                chain.add(NMEres)
                NMEres = None

    # 残基の順序を変更
    for model in structure:
        for chain in model:
            # 残基をリストとして取得してソート（番号＋挿入コード）
            residues = list(chain)
            residues.sort(key=lambda res: (res.id[1], res.id[2]))

            # チェーンからすべての残基を削除
            for res in list(chain):
                chain.detach_child(res.id)

            # ソート済み残基を再追加
            for res in residues:
                chain.add(res)


    # 変更後のPDBファイルを書き出す
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(output_dir)