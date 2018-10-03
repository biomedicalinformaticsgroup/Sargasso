class a(object):
    @staticmethod
    def stat():
        print("stat")

    @staticmethod
    def stat2():
        a.stat()
        print("stat2")

    def nonstat(self):
        print("non_stat")
        self.stat()




# a().stat()
# a().nonstat()
a().stat2()