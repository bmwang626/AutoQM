import os
import tarfile


def read_cosmo_tab_result_from_tar(f):
    """
    Modified from Yunsie's code
    """
    each_data_list = []
    # initialize everything
    solvent_name, solute_name, temp = None, None, None
    result_values = None
    line = f.readline()
    while line:
        # get the temperature and mole fraction
        if b"Settings  job" in line:
            temp = (
                line.split(b"T=")[1].split(b"K")[0].strip().decode("utf-8")
            )  # temp in K

        # get the result values
        if b"Nr Compound" in line:
            line = f.readline()
            solvent_name = line.split()[1].decode("utf-8")
            line = f.readline()
            solute_name = line.split()[1].decode("utf-8")
            result_values = line.split()[
                2:6
            ]  # H (in bar), ln(gamma), pv (vapor pressure in bar), Gsolv (kcal/mol)
            result_values = [
                result_value.decode("utf-8") for result_value in result_values
            ]
            # save the result as one list
            each_data_list.append(
                [solvent_name, None, solute_name, None, temp] + result_values + [None]
            )
            # initialize everything
            solvent_name, solute_name, temp = None, None, None
            result_values = None
        line = f.readline()
    return each_data_list


def get_dHsolv_value(each_data_list):
    # compute solvation enthalpy
    dGsolv_temp_dict = {}
    ind_298 = None
    for z in range(len(each_data_list)):
        temp = each_data_list[z][4]
        dGsolv = each_data_list[z][8]
        dGsolv_temp_dict[temp] = dGsolv
        if temp == "298.15":
            ind_298 = z
    if (
        "298.15" in dGsolv_temp_dict.keys()
        and "299.15" in dGsolv_temp_dict.keys()
        and "297.15" in dGsolv_temp_dict.keys()
    ):
        dGsolv_298 = float(dGsolv_temp_dict["298.15"])
        dSsolv_298 = -(
            float(dGsolv_temp_dict["299.15"]) - float(dGsolv_temp_dict["297.15"])
        ) / (299.15 - 297.15)
        dHsolv_298 = dGsolv_298 + 298.15 * dSsolv_298
        each_data_list[ind_298][9] = "%.8f" % dHsolv_298
    return each_data_list


def cosmo_parser(tar_file_path):
    each_data_lists = []
    try:
        tar = tarfile.open(tar_file_path)
    except tarfile.ReadError:
        print("tar file open failed")
        print(tar_file_path)
        return None
    try:
        for member in tar:
            if member.name.endswith(".tab"):
                f = tar.extractfile(member)
                each_data_list = read_cosmo_tab_result_from_tar(f)
                try:
                    each_data_list = get_dHsolv_value(each_data_list)
                except:
                    print("dHsolv calculation failed")
                    print(tar_file_path)
                    print(each_data_list)
                each_data_lists.append(each_data_list)
        tar.close()
    except tarfile.ReadError:
        print("tar file read failed")
        print(tar_file_path)
        return None
    return each_data_lists
