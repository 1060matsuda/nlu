# %%
import yaml
with open("./config.yaml", "r") as yml:
    config = yaml.safe_load(yml)

f1 = config["f1"]
"""GHz"""
f2 = config["f2"]
"""GHz"""
timestep = config["timestep"]

base_file="in.lammps.nonrb.shg"
modified_file="./mod_"+base_file
# %%
script = "#!/bin/bash\n"
script=script+"cp $MYWORK/common_files/"+base_file+" "+modified_file+"\n"
script=script+"sed -i -e 's/TIMESTEP_REPLACE_YAML/" + str(timestep)+"/g' "+modified_file+"\n"
script=script+"sed -i -e 's/F1_REPLACE_YAML/" + str(f1)+"/g' "+modified_file+"\n"
script=script+"sed -i -e 's/F2_REPLACE_YAML/" + str(f2)+"/g' "+modified_file

# %%
shell_script_file_name = "lammps_script_modifier_temp.sh"
"""with open(shell_script_file_name, mode="w") as f:
    f.write(script)"""
print(script)
#script should be catted and pasted in **.sh shell script.
