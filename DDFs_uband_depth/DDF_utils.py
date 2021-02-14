#Slight modification from the ddfInfo script in the opsimUtils.py module from the LSST_OpSim respository from Gordon Ricahrd's group. The change is that the function allows for the return of multiple proposal IDs matching the DDFs name. The most important part of this is that the EDFS is split into two subfields in FBS 1.7, and this is the only way to identify both of them simultaneously. 

# DDF RA/DEC dict
ddfCoord = {
    'COSMOS': (150.11, 2.14),
    'ELAISS1': (9.487, -44.0),
    'XMM-LSS': (35.707, -4.72),
    'ECDFS': (53.15, -28.08),
    '290': (349.377, -63.32),
    'EDFS': (61.28, -48.42)
}

def my_ddfInfo(opsimdb, ddfName):
    """
    Return DDF metainfo given the name and a opsim database object.

    Args:
        opsimdb: An opsim database object.
        ddfName(str): The name of the requested DDF field, e.g., COSMOS

    Returns:
        ddfInfo(dict): A dictionary containing metainfo (proposalId, RA/DEC and etc.) 
            for the requested DDF field. 
    """

    ddfName = str(ddfName)

    if ddfName not in ddfCoord.keys():
        print('DDF name provided is not correct! Please use one of the below: \n')
        print(list(ddfCoord.keys()))
        return None
    elif len(opsimdb.fetchPropInfo()[1]['DD']) == 0:
        print('No DDF in this Opsim run!')
        return None
    else:
        ddfInfo = {}
        propInfo = opsimdb.fetchPropInfo()[0]
        ddfInfo['proposalId'] = [key for key, elem in propInfo.items() if
                                 elem[3:3+len(ddfName)] == '{}'.format(ddfName)]
        ddfInfo['Coord'] = ddfCoord[ddfName]
        return ddfInfo
