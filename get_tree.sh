#!/bin/bash
(git ls-files && git ls-files --others --exclude-standard) | tree --fromfile
