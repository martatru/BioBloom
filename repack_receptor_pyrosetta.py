from pyrosetta import init, pose_from_pdb, get_fa_scorefxn
from pyrosetta.rosetta.core.pack.task import TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import RestrictToRepacking, InitializeFromCommandline, IncludeCurrent
from pyrosetta.rosetta.protocols.minimization_packing import PackRotamersMover

# uruchomienie PyRosetta z rozszerzonymi rotamerami Dunbrack
init("-mute all -ex1 -ex2aro -use_input_sc")

# wczytaj oczyszczony ACE
pose = pose_from_pdb("/home/marta/Desktop/docking_files/ACE_structures/ACE_clean.pdb")

# standardowa funkcja scoringu
scorefxn = get_fa_scorefxn()

# definiujemy zadanie: repacking tylko bocznych łańcuchów
tf = TaskFactory()
tf.push_back(InitializeFromCommandline())
tf.push_back(IncludeCurrent())        # uwzględnij oryginalne konformacje
tf.push_back(RestrictToRepacking())   # nie mutujemy, tylko zmieniamy rotamery

# generujemy i stosujemy mover
task = tf.create_task_and_apply_taskoperations(pose)
packer = PackRotamersMover(scorefxn, task)
packer.apply(pose)

# zapisujemy zrekonstruowany receptor
pose.dump_pdb("/home/marta/Desktop/docking_files/ACE_structures/ACE_repacked.pdb")
print("✅ Zapisano /home/marta/Desktop/docking_files/ACE_structures/ACE_repacked.pdb")
