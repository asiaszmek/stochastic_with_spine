from lxml import etree
#import argparse
import sys



def read_in_files(flist):
    roots = []
    species = set()
    specie_diff = {}
    rxn_ids = set()
    my_rxn_file = etree.Element("ReactionScheme")

    for fname_xml in flist:
        print(fname_xml)
        tree = etree.parse(fname_xml)
        roots.append(tree.getroot())
        for son in roots[-1]:
            if son.tag == "Specie":
                specie_id = son.attrib["id"]
                specie_kdiff = son.attrib["kdiff"]
                species.add(specie_id)
                if specie_id in specie_diff:
                    if specie_diff[specie_id]!= specie_kdiff:
                        print(specie_id, specie_diff[specie_id], specie_kdiff)
                else:
                    specie_diff[specie_id] = specie_kdiff
                    my_rxn_file.append(son)
            elif son.tag == "Reaction":
                if son.attrib["id"] in rxn_ids:
                    print(son.attrib["id"], " already exists", fname_xml)
                else:
                    rxn_ids.add(son.attrib["id"])
                    my_rxn_file.append(son)
    return my_rxn_file

flist_no_ER_spine_basal = [
    "Rxn_module_Ca.xml",
    "Rxn_module_RyR2CaM_KeizerSmith.xml",
    "Rxn_module_SERCA2.xml",
    "Rxn_module_SOCE.xml",
    ]

flist_no_ER_spine_basal_old_age = [
    "Rxn_module_Ca_old_age.xml",
    "Rxn_module_RyR2CaM_KeizerSmith.xml",
    "Rxn_module_SERCA2.xml",
    "Rxn_module_RyR2_KeizerSmith.xml",
    "Rxn_module_SOCE.xml",
    ]

flist_ER_spine_basal = flist_no_ER_spine_basal + [
    "Rxn_module_RyR3CaM_KeizerSmith.xml",
    "Rxn_module_SERCA3.xml",
    ]

flist_ER_spine_basal_old_age = flist_no_ER_spine_basal_old_age + [
    "Rxn_module_RyR3CaM_KeizerSmith.xml",
    "Rxn_module_SERCA3.xml",
    "Rxn_module_RyR3.xml",
    ]
flist_ER_spine_stim = flist_ER_spine_basal + ["Rxn_module_CaCbuf.xml"]
    
flist_no_ER_spine_stim = flist_no_ER_spine_basal + ["Rxn_module_CaCbuf.xml"]

flist_ER_spine_stim_old_age = flist_ER_spine_basal_old_age + ["Rxn_module_CaCbuf.xml"]
    
flist_no_ER_spine_stim = flist_no_ER_spine_basal_old_age + ["Rxn_module_CaCbuf.xml"]


flist_ER_spine_stim_Fluo4FF = flist_ER_spine_basal + ["Rxn_module_CaCbuf.xml",
    "Rxn_module_Fluo4FF.xml",
    ]
flist_no_ER_spine_stim_Fluo4FF = flist_no_ER_spine_basal + [
    "Rxn_module_CaCbuf.xml", "Rxn_module_Fluo4FF.xml"]


flist_ER_spine_stim_Fluo4FF_old_age = flist_ER_spine_basal_old_age + ["Rxn_module_CaCbuf.xml",
    "Rxn_module_Fluo4FF.xml",
    ]
flist_no_ER_spine_stim_Fluo4FF_old_age = flist_no_ER_spine_basal_old_age + [
    "Rxn_module_CaCbuf.xml", "Rxn_module_Fluo4FF.xml"]

if __name__ == "__main__":
    
    # 1 no mGluR no RyR
   

    my_rxn_f = read_in_files(flist_no_ER_spine_basal)
    with  open("Rxn_no_spine_ER_bas.xml", "w") as f:
        f.write(etree.tostring(my_rxn_f, pretty_print=True).decode("utf-8"))
    my_rxn_f = read_in_files(flist_ER_spine_basal)
    with  open("Rxn_spine_ER_bas.xml", "w") as f:
        f.write(etree.tostring(my_rxn_f, pretty_print=True).decode("utf-8"))
    my_rxn_f = read_in_files(flist_no_ER_spine_stim)
    with  open("Rxn_no_spine_ER.xml", "w") as f:
        f.write(etree.tostring(my_rxn_f, pretty_print=True).decode("utf-8"))
    my_rxn_f = read_in_files(flist_ER_spine_stim)
    with  open("Rxn_spine_ER.xml", "w") as f:
        f.write(etree.tostring(my_rxn_f, pretty_print=True).decode("utf-8"))
    my_rxn_f = read_in_files(flist_no_ER_spine_stim_Fluo4FF)
    with  open("Rxn_no_spine_ER_Fluo4FF.xml", "w") as f:
        f.write(etree.tostring(my_rxn_f, pretty_print=True).decode("utf-8"))
    my_rxn_f = read_in_files(flist_ER_spine_stim_Fluo4FF)
    with  open("Rxn_spine_ER_FLuo4FF.xml", "w") as f:
        f.write(etree.tostring(my_rxn_f, pretty_print=True).decode("utf-8"))


    my_rxn_f = read_in_files(flist_no_ER_spine_basal_old_age)
    with  open("Rxn_no_spine_ER_bas_old_age.xml", "w") as f:
        f.write(etree.tostring(my_rxn_f, pretty_print=True).decode("utf-8"))
    my_rxn_f = read_in_files(flist_ER_spine_basal)
    with  open("Rxn_spine_ER_bas_old_age.xml", "w") as f:
        f.write(etree.tostring(my_rxn_f, pretty_print=True).decode("utf-8"))
    my_rxn_f = read_in_files(flist_no_ER_spine_stim)
    with  open("Rxn_no_spine_ER_old_age.xml", "w") as f:
        f.write(etree.tostring(my_rxn_f, pretty_print=True).decode("utf-8"))
    my_rxn_f = read_in_files(flist_ER_spine_stim)
    with  open("Rxn_spine_ER_old_age.xml", "w") as f:
        f.write(etree.tostring(my_rxn_f, pretty_print=True).decode("utf-8"))
    my_rxn_f = read_in_files(flist_no_ER_spine_stim_Fluo4FF)
    with  open("Rxn_no_spine_ER_Fluo4FF_old_age.xml", "w") as f:
        f.write(etree.tostring(my_rxn_f, pretty_print=True).decode("utf-8"))
    my_rxn_f = read_in_files(flist_ER_spine_stim_Fluo4FF)
    with  open("Rxn_spine_ER_FLuo4FF_old_age.xml", "w") as f:
        f.write(etree.tostring(my_rxn_f, pretty_print=True).decode("utf-8"))
