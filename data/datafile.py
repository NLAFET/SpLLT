import vextractor
import os

class Datafile(object):

    def __init__(self, filename):

        # set filename
        self.filename = filename

        # open file if it exists
        if os.path.exists(filename):
            self.f = open(filename)
        else:
            # FIXME: handle error 
            print "File does not exist!"

    def __del__(self):

        # close file
        # FIXME: if file exists!
        self.f.close()
            
    # in datafile, get value corresponding to regular exp
    def get_value(self, exp):
        v = vextractor.get_value(self.f, exp)
        # position pointer to begining
        self.f.seek(0)

        return v
