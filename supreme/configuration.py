import yaml

from .utils import op_code_to_str


class Configuration(object):
    """
    This is a stub that just loads a yaml and turns itself into a dict-like object.

    In the future, this will add more features on checking fields, etc.
    """
    def __init__(self, configfile):
        """
        """

        with open(configfile) as f:
            self._config = yaml.load(f, Loader=yaml.SafeLoader)

    def patch_input_filename(self, filter_name, tract, patch_name):
        return '%s_%05d_%s_%s_patch_inputs.hs' % (self.outbase,
                                                  tract, patch_name, filter_name)

    def patch_map_filename(self, filter_name, tract, patch_name, map_type, operation):
        return "%s_%05d_%s_%s_%s_%s.hs" % (self.outbase,
                                           tract,
                                           patch_name,
                                           filter_name,
                                           map_type,
                                           op_code_to_str(operation))

    def tract_map_filename(self, filter_name, tract, map_type, operation):
        return "%s_%05d_%s_%s_%s.hs" % (self.outbase,
                                        tract,
                                        filter_name,
                                        map_type,
                                        op_code_to_str(operation))

    def __getattr__(self, attr):
        try:
            return self._config[attr]
        except KeyError:
            return object.__getattribute__(self, attr)

    def __setitem__(self, key, item):
        self._config.__setitem__(key, item)

    def __getitem__(self, key, item):
        return self._config.__getitem__(key, item)

    def __repr__(self):
        return self._config.__repr__()

    def __len__(self):
        return self._config.__len__()

    def __delitem__(self, key):
        self._config.__delitem__(key)

    def update(self, *args, **kwargs):
        return self._config.update(*args, **kwargs)

    def keys(self):
        return self._config.keys()

    def values(self):
        return self._config.values()

    def items(self):
        return self._config.items()

    def pop(self, key):
        return self._config.pop(key)

    def __cmp__(self, dict_):
        return self._config.__cmp__(dict_)

    def __contains__(self, item):
        return self._config.__contains__(item)

    def __iter__(self):
        return self._config.__iter__()

    def __unicode__(self):
        return self._config.__unicode__()
