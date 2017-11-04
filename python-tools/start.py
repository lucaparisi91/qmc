import input

# create the directories to study with changing number of walkers
def walker_study(walkers_to_study):
    
    #input.run_dir='/home/lucap/storage/cqmc/test'
    input.default_files.append("save.xml")
    input.default_files.append("walkers.xml")
    for m in walkers_to_study:
        a=input.new()
        a.root.find("algorithm").find("calculation").find("mean_walkers").text=str(m)
        a.root.find("algorithm").find("calculation").find("delta_walkers").text=str(m/4)
        a.make_jastrows()
        a.save()

# requires a certain time step
def time_step_study(times_to_study):
    input.default_files.append("save.xml")
    input.default_files.append("walkers.xml")
    for m in times_to_study:
        a=input.new()
        a.root.find("algorithm").find("calculation").find("delta_tau").text=str(m)
        a.make_jastrows()
        a.save()
