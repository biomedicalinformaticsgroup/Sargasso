class Manager(object):
    @classmethod
    def _create(cls, *args):
        raise NotImplementedError('Need to implement in subclass')

    @classmethod
    def get(cls, *args):
        raise NotImplementedError('Need to implement in subclass')
