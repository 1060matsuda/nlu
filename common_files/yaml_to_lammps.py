# %%
import yaml
with open("./config.yaml", "r") as yml:
    config = yaml.safe_load(yml)

suffixes = {"freqmix": ".mix", "shg": ".shg"}

are_multi_detectors_arranged = False
if "multi_detec" in config:
    # multi detector setting applied
    if config["multi_detec"] == True:
        are_multi_detectors_arranged = True

if "nlu_method" in config:
    nlu_method = str(config["nlu_method"])
else:
    print("WARNING: NLU method (freqmix / nlu) is not specified.")
    nlu_method = "shg"

suffix = suffixes[nlu_method]

f1 = config["f1"]
"""[GHz]"""
f2 = config["f2"]
"""[GHz]"""
timestep = config["timestep"]
"""[ps]"""
cycles_low = config["cycles_low"]
"""Input cycles of the lower frequency wave"""
# %%
if are_multi_detectors_arranged == False:
    tmp_script_name = "in.lammps"
else:
    tmp_script_name = "in.lammps.detec"
script_name = tmp_script_name+suffix
mod_script_name = "mod_"+script_name

script = "#!/bin/bash\n"
script = script+"cp $MYWORK/common_files/"+script_name+" ./"+mod_script_name+"\n"
script = script+"sed -i -e 's/TIMESTEP_REPLACE_YAML/" + \
    str(timestep)+"/g' ./"+mod_script_name+"\n"
script = script+"sed -i -e 's/CYCLES_LOW_REPLACE_YAML/" + \
    str(cycles_low)+"/g' ./"+mod_script_name+"\n"
script = script+"sed -i -e 's/F1_REPLACE_YAML/" + \
    str(f1)+"/g' ./"+mod_script_name+"\n"
script = script+"sed -i -e 's/F2_REPLACE_YAML/" + \
    str(f2)+"/g' ./"+mod_script_name



# %%
shell_script_file_name = "lammps_script_modifier_temp.sh"
"""with open(shell_script_file_name, mode="w") as f:
    f.write(script)"""
print(script)
# script should be catted and pasted in **.sh shell script.
