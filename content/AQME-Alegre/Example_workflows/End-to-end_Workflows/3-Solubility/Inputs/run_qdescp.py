from aqme.qdescp import qdescp

sdf_rdkit_files = f'CSEARCH/*.sdf'

qdescp(files=sdf_rdkit_files, boltz=True, program='xtb')