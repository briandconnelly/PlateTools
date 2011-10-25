# -*- coding: utf-8 -*-
# TODO: documentation

__author__ = "Brian Connelly <bdc@msu.edu>"
__version__ = 0.1

import numpy
from scipy import stats
import re

from PlateTools.Read import *
from PlateTools.Plate import *

class SMPNote(object):
    def __init__(self, note):
        self.note = note
    def __str__(self):
        return self.note

class SMPGroupSample(object):
    def __init__(self, name, group):
        self.name = name
        self.group = group
        self.values = []
    def __str__(self):
        num_values = len(self.values)
        valstring = "value"
        if num_values != 1:
            valstring += 's'
        return "Group '{gn}' Sample '{sn}' with {v} {vs}".format(gn=self.group.name, sn=self.name, v=num_values, vs=valstring)
    def add_value(self, val_tuple):
        self.values.append(val_tuple)
    def print_information(self):
        allvals = []

        print(self)
        for (well, value) in self.values:
            allvals.append(value)
            print("\tWell: {well}\tValue: {value}".format(well=well, value=value))
        print("")
        print("\tMean: {m}".format(m=numpy.mean(allvals)))
        print("\tStandard Deviation: {sd}".format(sd=numpy.std(allvals)))
        print("\tStandard Error: {se}".format(se=stats.sem(allvals)))
        print("")

class SMPGroup(object):
    def __init__(self, name):
        self.name = name
        self.samples = []
    def __str__(self):
        num_samples = len(self.samples)
        samplestring = "sample"
        if num_samples != 1:
            samplestring += 's'
        return "Group '{g}' with {s} {ss}".format(g=self.name, s=num_samples, ss=samplestring)
    def add_sample(self, sample):
        self.samples.append(sample)
    def print_information(self):
        print(self)
        for s in self.samples:
            num_values = len(s.values)
            valstring = "value"
            if num_values != 1:
                valstring += 's'
            print("\tSample '{sn}' with {v} {vs}:".format(sn=s.name, v=num_values, vs=valstring))

            allvals = []
            for (well, value) in s.values:
                allvals.append(value)
                print("\t\tWell: {well}\tValue: {val}".format(well=well, val=value))
            print("")
            print("\t\tMean: {m}".format(m=numpy.mean(allvals)))
            print("\t\tStandard Deviation: {sd}".format(sd=numpy.std(allvals)))
            print("\t\tStandard Error: {se}".format(se=stats.sem(allvals)))
            print("")


class SMPPlateRead(Read):
    def __init__(self, plate):
        super(SMPPlateRead, self).__init__(plate)
        self.info['time'] = None
        self.info['temperature'] = None
        self.data = None

    def __str__(self):
        return "Read of plate '{0}': Time={1}, Temperature={2}".format(self.plate.name,
                                                                       self.info['time'],
                                                                       self.info['temperature'])

class SMPPlate(Plate):
    info_fields = {'name': 'Name',
                   'export_version': 'Export Version',
                   'export_format': 'Export Format',
                   'read_type': 'Read Type',
                   'data_mode': 'Data Mode',
                   'data_type': 'Data Type',
                   'pre_read_included': 'Pre Read Included',
                   'kinetic_points': 'Kinetic Points',
                   'time_pattern': 'Time Patern',
                   'interval_density': 'Interval Density',
                   'start_wavelength': 'Start Wavelength',
                   'end_wavelength': 'End Wavelength',
                   'wavelength_step': 'Wavelength Step',
                   'num_wavelengths': 'Number of Wavelengths',
                   'read_wavelengths': 'Read Wavelengths',
                   'first_column_read': 'First Column Read',
                   'num_columns_read': 'Number of Columns Read',
                   'num_wells': 'Number of Wells',
                   'excitation_wavelengths': 'Excitation Wavelengths',
                   'cutoff_option': 'Cutoff Option',
                   'cutoff_filters': 'Cutoff Filters',
                   'sweep_wave_option': 'Sweep Wave Option',
                   'sweep_fixed_wave': 'Sweep Fixed Wave',
                   'reads_per_well': 'Reads Per Well',
                   'pmt_setting': 'PMT Setting',
                   'start_integration_time': 'Start Integration Time',
                   'end_integration_time': 'End Integration Time',
                   'first_row_read': 'First Row Read',
                   'num_rows_read': 'Num Rows Read',
                   'time_tags': 'Time Tags'}

    def __init__(self, num_wells):
        super(SMPPlate, self).__init__(num_wells=num_wells)

    def __str__(self):
        num_reads = len(self.reads)
        s = "'{0}': {1} data with {2} {3} read".format(self.info['name'],
                                                       self.info['data_mode'],
                                                       num_reads,
                                                       self.info['read_type'])
        if num_reads != 1: s += 's'
        return s

    def print_information(self):
        print(self)
        for k, v in self.info_fields.iteritems():
            try:
                print("\t{desc}: {val}".format(desc=v, val=self.info[k]))
            except KeyError as err:
                print("\t{desc}: Unknown".format(desc=v))
        print("")

class Cuvette(object):
    def __init__(self):
        # TODO: implement
        pass

class Group(object):
    def __init__(self):
        pass

class Experiment(object):
    def __init__(self, fp):
        self.notes = []
        self.plates = {}
        self.cuvettes = {}
        self.groups = {}
        self.read_file(fp)

    def __str__(self):
        return "Experiment [{0} Notes], {1} Plates, {2} Cuvettes".format(len(self.notes), len(self.plates), len(self.cuvettes))

    def plate_names(self):
        return [p.info['name'] for p in self.plates]

    def cuvette_names(self):
        return [p.name for p in self.cuvettes]

    def read_block(self, fp):
        block_lines = []

        while True:
            last_pos = fp.tell()
            line = fp.readline()
            if not line: break

            line = line.rstrip()

            m = re.match('^(.*)~End', line)
            if m:
                if len(m.group(1)) >= 1:
                    block_lines.append(m.group(1))
                break
            else:
                block_lines.append(line)

        return block_lines

    def parse_note_block(self, block):
        """ Parse a Note block """
        note = ''

        for line in block[1:]:
            note += line.rstrip()

        N = SMPNote(note)
        return N

    def parse_group_block(self, block):
        """ Parse a Group block """
        for line in block:
            line_tokens = line.split('\t')

            if line_tokens[0] == "Group:":
                g = SMPGroup(line_tokens[1])
            elif line_tokens[0] == "Sample":
                s = None
            else:
                if len(line_tokens) == 5:
                    if s:
                        g.add_sample(s)
                    s = SMPGroupSample(name=line_tokens[0], group=g)
                    s.add_value((line_tokens[1], float(line_tokens[3])))
                elif len(line_tokens) == 4:
                    s.add_value((line_tokens[1], float(line_tokens[3])))
                elif line == '':
                    g.add_sample(s)
                else:
                    print("ERROR: Unexpected line")

        return g

    def parse_cuvette_block(self, block):
        """ Parse a Cuvette block """
        # TODO: implement this.  We just have header infos now.
        for line in block:
            m = re.match("^Cuvette:", line)
            if m:
                information = line.split('\t')
                name = information[1]
                export_version = information[2]
                export_format = information[3]

                read_type = information[4]
                if read_type not in ['Endpoint', 'Kinetic', 'Spectrum', 'Well Scan', 'Flex']:
                    print("ERROR: Invalid data format (read_type)")

                data_mode = information[5]
                data_type = information[6]
                pre_read_included = information[7] == "TRUE"
                kinetic_points = int(information[8])
                time_pattern = information[9]
                interval_density = information[10]
                start_wavelength = information[11]
                if start_wavelength: start_wavelength = float(start_wavelength)
                end_wavelength = information[12]
                if end_wavelength: end_wavelength = float(end_wavelength)
                wavelength_step = information[13]
                if wavelength_step: wavelength_step = float(wavelength_step)
                num_wavelengths = information[14]
                read_wavelengths = information[15] # TODO: int??
                first_column_read = int(information[16])
                num_columns_read = int(information[17])
                num_cuvettes = int(information[18])

                if data_mode == "Fluorescence" or data_mode == "Luminescence":
                    excitation_wavelengths = information[19] # TODO: space separated.  Parse into list. of floats?

                    cutoff_option = information[20]
                    if cutoff_option not in ['Manual', 'Automatic', '']:
                        print("ERROR: Invalid data format (cutoff_option)")

                    cutoff_filters = information[21] # TODO: space separated.  Parse into list. of floats??

                    sweep_wave_option = information[22] # TODO validate
                    if sweep_wave_option not in ['EmSweep', 'ExSweep', '']:
                        print("ERROR: Invalid data format (sweep_wave_option)")

                    sweep_fixed_wave = information[23]

                continue

    def parse_plate_block(self, block):
        """ Parse a Plate block """
        info_read = False
        header_read = False
        read_data = False

        for line in block:
            m = re.match("^Plate:", line)
            if m:
                info = {}
                info_raw = line.split('\t')

                info['name'] = info_raw[1]
                info['export_version'] = info_raw[2]
                info['export_format'] = info_raw[3]

                info['read_type'] = info_raw[4]
                if info['read_type'] not in ['Endpoint', 'Kinetic', 'Spectrum', 'Well Scan', 'Flex']:
                    print("ERROR: Invalid data format (read_type)")

                info['data_mode'] = info_raw[5]
                info['data_type'] = info_raw[6]
                info['pre_read_included'] = info_raw[7] == "TRUE"
                info['kinetic_points'] = int(info_raw[8])
                info['time_pattern'] = info_raw[9]
                info['interval_density'] = info_raw[10]
                info['start_wavelength'] = info_raw[11]
                if info['start_wavelength']: info['start_wavelength'] = float(info['start_wavelength'])
                info['end_wavelength'] = info_raw[12]
                if info['end_wavelength']: info['end_wavelength'] = float(info['end_wavelength'])
                info['wavelength_step'] = info_raw[13]
                if info['wavelength_step']: info['wavelength_step'] = float(info['wavelength_step'])
                info['num_wavelengths'] = info_raw[14]
                info['read_wavelengths'] = info_raw[15] # TODO: int??
                info['first_column_read'] = int(info_raw[16])
                info['num_columns_read'] = int(info_raw[17])
                info['num_wells'] = int(info_raw[18])

                if info['data_mode'] == "Fluorescence" or info['data_mode'] == "Luminescence":
                    info['excitation_wavelengths'] = info_raw[19] # TODO: space separated.  Parse into list. of floats?

                    info['cutoff_option'] = info_raw[20]
                    if info['cutoff_option'] not in ['Manual', 'Automatic', '']:
                        print("ERROR: Invalid data format (cutoff_option)")

                    info['cutoff_filters'] = info_raw[21] # TODO: space separated.  Parse into list. of floats??

                    info['sweep_wave_option'] = info_raw[22] # TODO validate
                    if info['sweep_wave_option'] not in ['EmSweep', 'ExSweep', '']:
                        print("ERROR: Invalid data format (sweep_wave_option)")

                    info['sweep_fixed_wave'] = info_raw[23]
                    info['reads_per_well'] = int(info_raw[24])

                    info['pmt_setting'] = info_raw[25] # TODO: validate
                    if info['pmt_setting'] not in ['Automatic', 'High', 'Medium', 'Low']:
                        print("ERROR: Invalid data format (pmt_setting)")

                    info['start_integration_time'] = info_raw[26]
                    info['end_integration_time'] = info_raw[27] # TODO: float??
                    info['first_row_read'] = int(info_raw[28])
                    info['num_rows_read'] = int(info_raw[29])
                    info['time_tags'] = info_raw[30] # TODO: boolean?? - dunno, but value can be 'None'
                else:
                    info['first_row_read'] = int(info_raw[19])
                    info['num_rows_read'] = int(info_raw[20])
                    info['time_tags'] = info_raw[21] # TODO:boolean?? - dunno, but value can be 'None'

                info_read = True

                P = SMPPlate(num_wells=info['num_wells'])
                P.info = info
                P.info['raw'] = info_raw

                continue

            m = re.match("^(Time\(hh:mm:ss\))?\tTemperature", line)
            if m and info_read:
                # Get the number of rows and columns in the plate.  These
                # should usually be able to calculated assuming 2:3, but just
                # in case...
                P.info['num_columns_plate'] = len(line.split('\t')) - 2 # -2 because first two are time/temp
                P.info['num_rows_plate'] = P.info['num_wells'] / P.info['num_columns_plate'] # Is this a valid computation?

                P.num_rows = P.info['num_rows_plate']
                P.num_columns = P.info['num_columns_plate']
                header_read = True
                continue

            if info_read and header_read:
                line_tokens = line.split('\t')

                if len(line_tokens) >= 2 and (line_tokens[0] or line_tokens[1]):
                    R = SMPPlateRead(P)
                    data = []
                    data_row_num = 1
                    bad_data = False
                else:
                    data_row_num += 1

                if len(line_tokens) >= 1 and line_tokens[0]:
                    R.info['time'] = line_tokens[0]

                if len(line_tokens) >= 2 and line_tokens[1]:
                    R.info['temperature'] = float(line_tokens[1])

                if data_row_num >= P.info['first_row_read'] and data_row_num < (P.info['first_row_read'] + P.info['num_rows_read']):
                    row_data = line_tokens[2+(P.info['first_column_read']-1):]
                    row_data = row_data[:P.info['num_columns_read']]
                    try:
                        row_data = map(float, row_data)
                        data.append(row_data)
                    except ValueError:
                        bad_data = True
                elif data_row_num <= P.num_rows:
                    # This is a valid plate row, but it wasn't read
                    pass
                else:
                    # We've reached the end of a read.  Append the data to the set of reads if bad_data==False
                    if not bad_data:
                        R.data = numpy.array(data)
                        P.add_read(R)
                continue

            else:
                print("ERROR: Unexpected line")

        return P

    def read_file(self, fp):
        """ Read a data file """

        while True:
            last_pos = fp.tell()
            line = fp.readline()
            if not line: break

            m = re.match('^\s*##BLOCKS=\s*(\d+)\s*', line)
            if m:
                num_blocks = int(m.group(1))

                for i in range(num_blocks):
                    b = self.read_block(fp)
                    if len(b) == 0:
                        print("ERROR: bad block")
                        return

                    m = re.match('^(Note|Plate|Cuvette|Group)', b[0])
                    if m:
                        if m.group(1) == 'Note':
                            note = self.parse_note_block(b)
                            self.notes.append(note)
                        elif m.group(1) == 'Plate':
                            plate = self.parse_plate_block(b)
                            self.plates[plate.info['name']] = plate
                        elif m.group(1) == 'Cuvette':
                            cuvette = self.parse_cuvette_block(b)
                            self.cuvettes[cuvette] = cuvette
                        elif m.group(1) == 'Group':
                            group = self.parse_group_block(b)
                            self.groups[group.name] = group
                    else:
                        print("ERROR: Invalid block type")

                continue

            m = re.match('^Copyright.*$', line)
            if m:
                self.copyright = line.rstrip()
                continue

            print("ERROR: Unexpected line")
