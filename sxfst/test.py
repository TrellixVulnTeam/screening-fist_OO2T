import os
import random
from unittest import TestCase
from sxfst import PlateData, Cpd, Screen

root = '../nb/00-pilot/00.3/platereader/'
paths = random.choices([os.path.join(root, i) for i in os.listdir(root)],
                       k=4)
picklist = '../nb/00-pilot/00.3/echo/picklists/2022-03-12-00.3-picklist.csv'

class TestPlateData(TestCase):
    def testInit(self):
        for path_ in paths:
            data = PlateData(path_)
    def testMeta(self):
        for path_ in paths:
            data = PlateData(path_)
            meta = data.metadata
    def testDf(self):
        for path_ in paths:
            data = PlateData(path_)
            df = data.df
    def testTimestamp(self):
        for path_ in paths:
            data = PlateData(path_)
            _ = data.timestamp
    def testIter(self):
        for path_ in paths:
            data = PlateData(path_)
            for i in data:
                pass
    def testGetitem(self):
        for path_ in paths:
            data = PlateData(path_)
            data[0]
            data[3:10]
            data['A1']
            data[['A1', 'A2','A4']]
    def testRepr(self):
        for path_ in paths:
            data = PlateData(path_)
            #print(data)
            str(data)
    def testLen(self):
        for path_ in paths:
            data = PlateData(path_)
            _ = len(data)

class TestCpd(TestCase):
    def testInit(self):
        # csvs
        #cpd = Cpd(csv='../nb/00-pilot/00.3/echo/picklists/2022-03-12-00.3-picklist.csv',
        #          search='A1')
        pass


class TestScreen(TestCase):
    def testInit(self):
        scrn = Screen(picklist=picklist,
                      files=paths)
    def testLen(self):
        scrn = Screen(picklist=picklist,
                      files=paths)
        len(scrn)
    def testGetitem(self):
        scrn = Screen(picklist=picklist,
                      files=paths)
        #print(scrn['S4206'])
        scrn[9]
        scrn[-1]
        scrn[100:200]
