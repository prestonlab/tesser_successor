import pandas


# receives a tab-separated input file and returns the data frame containing its contents
def read_file (directory):
    with open (directory) as file:
        matrix = pandas.read_csv (file, sep='\t')
    return matrix

# gets the objects column from the experimental run which is inputted
def get_objects(directory):
    matrix = read_file (directory)
    obj_sequence = matrix[[' objnum']]    
    return obj_sequence
