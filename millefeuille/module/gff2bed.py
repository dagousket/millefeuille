import os
import sys
import argparse
import re


def get_featureDict(featureInfo: str) -> dict:
    """
    Create a dictionary from the feature field of a GFF file.

    The object "featureInfo" should be the last field of a tab-delimited GFF file.
    And it really must be a GFF3 file, not a GTF or some such nonsense.

    Parameters
    ----------
    featureInfo
      feature field of a GFF file, which contains the information of the features

    Returns
    -------
    dict
      A dictionary with the feature type as key and a list of values.
    """
    featureDict = {}
    atoms = [k.lstrip(" ").rstrip(" ") for k in featureInfo.split(";")]
    for k in atoms:
        index, entry = k.split("=")
        featureDict[index.rstrip(" ").lstrip(" ")] = [
            j.rstrip(" ").lstrip(" ") for j in entry.split(",")
        ]
    return featureDict


def get_Dictgff(
    file_gff: str,
    mol_type: str = "exon",
    feature_type: str = "Parent",
    id_as_features: bool = True,
) -> dict:
    """
    Create a dictionary from the GFF file which contains the infomations of all feature_type of a given mol_type (ie exon, CDS).

     Parameters
     ----------
     file_gff
       The GFF file to be converted
     mol_type
       The molecular type (column 3 of the GFF file) selected for the BED files, default is exon
     feature_type
       The feature type (column 9 of the GFF file) selected for the BED files, default is Parent
     id_as_features
       If set to True, the ID of each element will be set as a string containing all its features

     Returns
     -------
     dict
       A dictionary with the feature type as key and a dictionary for each mol_type in the feature with chr, start, stop and strand keys as value.
    """
    with open(file_gff, "r") as f:
        myDict = {}
        for l in f:
            line = l.rstrip("\r\n").split("\t")
            if len(line) > 3 and line[2] == mol_type:
                # get feature type in key of featDict
                featDict = get_featureDict(line[-1])
                # get mol_type info
                keys = ["chr", "start", "stop", "strand"]
                values = [line[i] for i in [0, 3, 4, 6]]
                # add the key name if id_as_feature argument is specified
                if id_as_features:
                    keys.append("name")
                    values.append(line[-1])
                exonDict = dict(zip(keys, values))
                # create dict key=feattype value=moltypeDict1
                for t in featDict[feature_type]:
                    if t in myDict:
                        # add the exon
                        myDict[t].append(exonDict)
                    else:
                        # create key and add exon
                        myDict[t] = [exonDict]
    f.closed
    return myDict


def consistency_check(
    myDict: dict,
    feature_type: str = "Parent",
    verbose: bool = True,
    mol_type: str = "exon",
    discard: bool = True,
) -> bool:
    """
    Check the consistency of the data from gff.
           Return false if chromosome or strand consistency is not respected.
           Return true if respected and order the list of mol_type dictionary according to start postition.

     Parameters
     ----------
     myDict
       The dictionary containing the information of the features
     feature_type
       The feature type (column 9 of the GFF file) selected for the BED files, default is Parent
     verbose
       If set to True, the function will print warnings for strand inconsistencies and overlapping elements.
     mol_type
       The molecular type (column 3 of the GFF file) selected for the BED files, default is exon
     discard
       If set to True, the function will discard elements that raise warnings in the consistency check.

     Returns
     -------
     bool
       True if the data is consistent, False otherwise.
    """
    alloverlap = []
    allstrand = []
    allchr = 0

    for t in myDict:
        # chromosome consistency
        if len(set([d["chr"] for d in myDict[t]])) != 1:
            allchr += 1
            print(
                "WARNING : Chromosome consistency is not satisfied for "
                + feature_type
                + " "
                + t
            )
        # strand consistency
        elif len(set([d["strand"] for d in myDict[t]])) != 1:
            allstrand.append(t)
        # order exon based on start values
        startkey = [int(k["start"]) for k in myDict[t]]
        exonsort = [x for (y, x) in sorted(zip(startkey, myDict[t]))]
        myDict[t] = exonsort
        # check overlapping
        initstop = -1
        overlap = False
        for exon in myDict[t]:
            if int(exon["start"]) < int(initstop):
                overlap = True
                # print str(exon['start'])+' < '+str(initstop)+' in '+str(t)
                # Commented line print the position where exon overlapping is detected
            initstop = exon["stop"]
        if overlap:
            alloverlap.append(t)
    # Count and print warnings
    if alloverlap and verbose:
        print(
            "\nWarning : overlapping "
            + mol_type
            + " in "
            + feature_type
            + " :\n"
            + ",".join(alloverlap)
        )
    if allstrand and verbose:
        print(
            "\nWarning : Strand consistency is not satisfied for "
            + feature_type
            + " :\n"
            + t
            + ",".join(allstrand)
        )
    print(
        "\nTotal of :\n\t"
        + str(len(allstrand))
        + " strand inconsistencies\n\t"
        + str(len(alloverlap))
        + " "
        + feature_type
        + " with overlapping "
        + mol_type
        + "\n\t"
        + str(allchr)
        + " chromosome inconsistencies"
    )
    # Discard inconsistent element
    if discard:
        nbdisc = len(myDict)
        reject = list(set(allstrand) | set(alloverlap))
        for k in reject:
            myDict.pop(k)
        nbdisc -= len(myDict)
        print(
            "\nDiscarded "
            + str(nbdisc)
            + " elements with strand inconsistencies and/or overlapping\n"
        )
    return True


def bed6_generator(
    bedname: str,
    file_gff: str,
    mol_type: str = "exon",
    feature_type: str = "Parent",
    path: str = "./",
    id_as_features: bool = True,
    skip_exon_number: bool = True,
) -> None:
    """
    Create a simple BED file in the working directory.

    Parameters
       ----------
       file_gff
           The GFF file to be converted
       mol_type
           The molecular type (column 3 of the GFF file) selected for the BED files, default is exon
       feature_type
           The feature type (column 9 of the GFF file) selected for the BED files, default is Parent
       path
           The location where BED files will be created, default is current working directory
       bedname
           The name of the BED files, default is the GFF file name
       id_as_features
           Will set the ID of each element as a string containing all its features
       skip_exon_number
           If set, the program will skip adding _# for exon number.

       Returns
       -------
       None, but creates a BED6 file in the specified path.
    """
    # uses os module
    myDict = get_Dictgff(file_gff, mol_type, feature_type)
    if consistency_check(myDict):
        # only allow to pursue script if check script runs correctly
        with open(path.rstrip("/") + "/" + bedname + ".bed6", "w") as f6:  # check
            for t in myDict:
                exnum = 0
                for exon in myDict[t]:
                    # add a number to exons
                    exnum += 1
                    if id_as_features:
                        name = exon["name"]
                    else:
                        if skip_exon_number is False:
                            name = "_".join([t, str(exnum)])
                        else:
                            name = t
                    f6.write(
                        "\t".join(
                            [
                                exon["chr"],
                                str(int(exon["start"]) - 1),
                                exon["stop"],
                                name,
                                ".",
                                exon["strand"] + "\n",
                            ]
                        )
                    )
        f6.closed


def bed12_generator(
    bedname: str,
    file_gff: str,
    mol_type: str = "exon",
    feature_type: str = "Parent",
    path: str = "./",
    check: bool = False,
    id_as_features: bool = True,
) -> None:
    """
    Create a BED12 file in the working directory.

     Parameters
       ----------
       file_gff
           The GFF file to be converted
       mol_type
           The molecular type (column 3 of the GFF file) selected for the BED files, default is exon
       feature_type
           The feature type (column 9 of the GFF file) selected for the BED files, default is Parent
       path
           The location where BED files will be created, default is current working directory
       bedname
           The name of the BED files, default is the GFF file name
       check
           If set to True, the function will check the consistency of the data before creating the BED12 file.
       id_as_features
           Will set the ID of each element as a string containing all its features

       Returns
       -------
       None, but creates a BED12 file in the specified path.
    """
    # uses os module
    myDict = get_Dictgff(file_gff, mol_type, feature_type)
    if not check:
        if consistency_check(myDict):
            check = True
    if check:
        # only allow to pursue script if check script runs correctly
        with open(path.rstrip("/") + "/" + bedname + ".bed12", "w") as f12:
            for t in myDict:
                blockSizes = [int(k["stop"]) - int(k["start"]) for k in myDict[t]]
                blockStart = [
                    int(k["start"]) - int(myDict[t][0]["start"]) for k in myDict[t]
                ]
                # Handling size 0 exons appearing in gff file
                blockStart = [
                    blockStart[k] for k in range(len(blockSizes)) if blockSizes[k] > 0
                ]
                blockSizes = [
                    blockSizes[k] for k in range(len(blockSizes)) if blockSizes[k] > 0
                ]
                blockCount = str(len(blockSizes))
                # Handling single size 0 exon in a transcript
                if not blockSizes:
                    continue
                # Settinf int as string
                blockStart = ",".join(str(e) for e in blockStart)
                blockSizes = ",".join(str(e) for e in blockSizes)
                # check if correct by :
                # blockStart[-1]+myDict[t][0]['start']+blockSizes[-1] == myDict[t][-1]['stop']
                chromStart = str(int(myDict[t][0]["start"]) - 1)
                chromEnd = myDict[t][-1]["stop"]
                if id_as_features:
                    # recreate the list of feature without the feat_type from argument
                    featD = get_featureDict(myDict[t][0]["name"])
                    featD["Name"] = [re.sub(":[0-9]$", "", featD["Name"][0])]
                    name = ";".join(
                        [
                            str(k) + "=" + ",".join(v)
                            for (k, v) in featD.items()
                            if k != "Parent"
                        ]
                        + ["Parent=" + str(t)]
                    )
                else:
                    name = t
                f12.write(
                    "\t".join(
                        [
                            myDict[t][0]["chr"],
                            chromStart,
                            chromEnd,
                            name,
                            "0",
                            myDict[t][0]["strand"],
                            chromStart,
                            chromEnd,
                            "0",
                            blockCount,
                            blockSizes,
                            blockStart + "\n",
                        ]
                    )
                )
        f12.closed


def gff2bed(
    gff_file: str,
    bed12: bool = True,
    no_bed6: bool = False,
    mol_type: str = "exon",
    feature_type: str = "Parent",
    id_as_features: bool = True,
    path: str = "./",
    name: str = "noinp",
    verbose: bool = True,
    discard: bool = True,
    skip_exon_number: bool = True,
) -> None:
    """
    Creates BED files from GFF file. In BED12, groups all the elements of a selected molecular type according to their feature type.

    Parameters
       ----------
       gff_file
           Name of the GFF file to be converted
       bed12
           Creates the corresponding BED12 file
       no_bed6
           Prevents from creating the corresponding BED6 file
       mol_type
           The molecular type (column 3 of the GFF file) selected for the BED files, default is exon
       feature_type
           The feature type (column 9 of the GFF file) selected for the BED files, default is Parent
       id_as_features
           Will set the ID of each element as a string containing all its features
       path
           The location where BED files will be created, default is current working directory
       name
           The name of the BED files, default is the GFF file name
       verbose
           Will outpout in stdout the command arguments and the name of each element raising a warning in consistency check
       discard
           Will discard the element raising a warning in strand consistency and overlapping check
       skip_exon_number
           If set, the program will skip adding _# for exon number.

       Returns
       -------
       None, but creates BED files in the specified path.
    """

    if name == "noinp":
        bedname = os.path.basename(gff_file).replace(".gff", "")
    else:
        bedname = name
    check = False
    if no_bed6:
        print("\nCreating BED6 with all " + mol_type + "s from the gff file")
        bed6_generator(gff_file, mol_type, feature_type, path, bedname)
        # Avoid redoing consistency check in bed12 if already done for bed6
        check = True
    if bed12:
        print(
            "\nCreating BED12 with all "
            + mol_type
            + "s from gff file, grouped according to their "
            + feature_type
        )
        bed12_generator(gff_file, mol_type, feature_type, path, bedname, check)
