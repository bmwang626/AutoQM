#!/usr/bin/env python3
# encoding: utf-8

"""
Copied from https://github.com/xiaoruiDong/easy_rmg_model/blob/master/easy_rmg_model/template_writer/base.py
"""

import os.path

from jinja2 import Environment, FileSystemLoader, Template


class BaseTemplateWriter(object):
    """
    This is a base class for all template writer.
    """

    default_settings = {"template_file": None, "save_path": "./base_template.txt"}
    default_template = ""

    def __init__(self, spec={}):
        for key, value in spec.items():
            setattr(self, key, value)
        for key, value in self.default_settings.items():
            if not hasattr(self, key):
                setattr(self, key, value)

    @property
    def template(self):
        if self.template_file:
            template_dir, template_file = os.path.split(self.template_path)
            env = Environment(loader=FileSystemLoader(template_dir), autoescape=True)
            template = env.get_template(template_file)
        else:
            template = Template(self.default_template, autoescape=True)
        return template

    @property
    def rendered_template(self):
        if not hasattr(self, "_rendered"):
            if self.template:
                self._rendered = self.template.render(self.to_dict())
        return self._rendered

    def save(self):
        with open(self.save_path, "w") as f:
            f.write(self.rendered_template)
        return True

    def to_dict(self):
        pass

    @classmethod
    def from_dict(cls, spec_dict):
        return cls(spec_dict)
