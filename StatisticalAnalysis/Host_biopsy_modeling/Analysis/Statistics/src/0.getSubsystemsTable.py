# Porthmeus
# 11.02.22

# small script to extract the subsystems of the model annotation into a table
import cobra as cb
import pandas as pd

mod = cb.io.read_sbml_model("../../../resources/models/colormore3D.xml")

subsystems = []
reactions = []
for grp in mod.groups:
    subsystems.extend([grp.name] * len(grp.members))
    reactions.extend([x.id for x in grp.members])

df = pd.DataFrame({"subsystem" : subsystems, "reaction" : reactions})

df.to_csv("../dat/subsystems.csv", index = False)
