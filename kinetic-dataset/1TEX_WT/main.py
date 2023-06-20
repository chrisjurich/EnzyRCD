from pymol import cmd


cmd.load('1tex.cif')

print(cmd.get_fastastr())

cmd.delete('all')

cmd.load('production.pdb')
print(cmd.get_fastastr())
