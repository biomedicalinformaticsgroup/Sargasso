class Manager(object):
    @classmethod
    def _create(cls, data_type):
        raise NotImplementedError('Need to implement in subclass')

    @classmethod
    def get(cls, data_type):
        raise NotImplementedError('Need to implement in subclass')
