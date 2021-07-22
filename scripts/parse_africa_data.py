import numpy
import parse_midas_data
import config

def parse_subject_sample_map(sample_metadata_map = {}):

    import config

    if len(sample_metadata_map)==0:
        sample_metadata_map = parse_sample_metadata_map()


    subject_sample_map = {}
    for sample_id in sample_metadata_map:
        subject_id, dummy, accession_id, country, continent, order = sample_metadata_map[sample_id]

        if subject_id not in subject_sample_map:
            subject_sample_map[subject_id] = {}

        if sample_id not in subject_sample_map[subject_id]:
            subject_sample_map[subject_id][sample_id] = set()

        subject_sample_map[subject_id][sample_id].add(accession_id)


    return subject_sample_map




###############################################################################
#
# Loads metadata for HMP samples
# Returns map from sample -> (subject_id, sample_id, accession_id, country, continent, temporal_order)
#
###############################################################################
def parse_sample_metadata_map():

    sample_metadata_map = {}

    # First load Ethiopia metadata
    #file = open(config.scripts_directory+"HMP_ids_order.txt","r")
    file = open(config.scripts_directory+"PRJNA504891_run_accessions.txt","r")
    file.readline() # header
    for line in file:
        items = line.split("\t")
        sample_id = items[0].strip()
        accession_id = items[1].strip()
        # one sample per subject
        subject_id = sample_id
        country = "Ethiopia"
        continent = "Africa"
        order = 1

        sample_metadata_map[sample_id] = (subject_id, sample_id, accession_id, country, continent, order)


    file.close()


    # load Madagascar metadata
    file = open(config.scripts_directory+"PRJNA485056_run_accessions.txt","r")
    file.readline() # header
    for line in file:
        items = line.split("\t")
        sample_id = items[0].strip()
        accession_id = items[1].strip()
        # one sample per subject
        subject_id = sample_id
        country = "Madagascar"
        continent = "Africa"
        order = 1

        sample_metadata_map[sample_id] = (subject_id, sample_id, accession_id, country, continent, order)

    file.close()

    return sample_metadata_map
