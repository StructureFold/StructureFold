#!/usr/bin/python

import threading


class subThread (threading.Thread):
    def __init__(self, threadID, function, *para):


        threading.Thread.__init__(self)
        self.para = para
        self.fun = function
        self.threadID = threadID
        self.name = "Thread-"+self.threadID
    def run(self):
        print "Starting " + self.name
        self.fun(*self.para)
        print "Exiting " + self.name


