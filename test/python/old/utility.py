from nose.tools import assert_almost_equals
import itertools

def assert_almost_equals_lists(list1, list2, places = 6):
     
    '''Assert that the contents of List1 and List2 are
    equal to the supplied number of decimal places. Do not support nested
    lists of depth more than 2.
    '''
     
    # If the list is nested open it out
    if type(list1[0]).__name__ == "list" or type(list1[0]).__name__ == "tuple":
         
        for entry1, entry2 in itertools.izip_longest(list1, list2):
            placesList = [places] * len(entry1)
            map(assert_almost_equals, entry1, entry2, placesList)     
    else:
        places = [places] * len(list1)
        map(assert_almost_equals, list1, list2, places)