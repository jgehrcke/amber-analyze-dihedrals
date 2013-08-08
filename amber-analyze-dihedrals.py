#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#   Copyright (C) 2013 Jan-Philip Gehrcke
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.


import os
import sys
import json
import yaml
import codecs
import argparse
import itertools
import operator
from subprocess import Popen, PIPE
import logging


logging.basicConfig(
    format='%(asctime)s,%(msecs)-6.1f %(levelname)s: %(message)s',
    datefmt='%H:%M:%S')
log = logging.getLogger()
log.setLevel(logging.INFO)


class Ambmask(object):
    """Ambmask object prepares ambmask runs with equal Amber topology and
    coordinate file and printlevel/outformat settings.
    """
    def __init__(
            self,
            topology_file_path,
            coordinate_file_path,
            printlevel="0",
            outformat="amber",
            ):
        assert os.path.isfile(topology_file_path)
        assert os.path.isfile(coordinate_file_path)
        self._topology_file_path = topology_file_path
        self._coordinate_file_path = coordinate_file_path
        self._printlevel = printlevel
        self._outformat = outformat

    def _run(self):
        """Run subprocess. If ambmask exits with return code other than 0,
        exit this program. Otherwise return ambmask's stdout.
        """
        args = ["ambmask",
            "-p",
            self._topology_file_path,
            "-c",
            self._coordinate_file_path,
            "-prnlev",
            self._printlevel,
            "-out",
            self._outformat,
            "-find",
            self._maskstring
            ]
        try:
            sp = Popen(args, stdout=PIPE, stderr=PIPE)
        except OSError as e:
            sys.exit(("Error executing ambpdb. Is it in your PATH? Specific "
                "error: '%s'" % e))
        out, err = sp.communicate()
        returncode = sp.returncode
        if returncode != 0:
            sys.stderr.write("ambmask error.\n")
            sys.stderr.write("ambmask stderr: %s\n" % err)
            sys.stderr.write("ambmask stdout: %s\n" % out)
            sys.stderr.flush()
            sys.exit(1)
        return out

    def residue_ids_by_name(self, resname):
        """For a given residue name, return a list of matching residue IDs.
        """
        self._maskstring = ":%s" % resname
        ambmask_stdout = self._run()
        if not ambmask_stdout:
            log.debug("%s did not write stdout." % self)
            return None
        residue_ids = []
        for rawline in ambmask_stdout.splitlines():
            line = rawline.strip()
            if line.startswith("RES"):
                # output:
                # RES        2       2
                # RES        4       4
                tokens = line.split()
                start = int(tokens[1])
                end = int(tokens [2])
                if start == end:
                    residue_ids.append(start)
                elif not end > start:
                    sys.exit("Expected %s to be larger than %s in line %s" %
                        (end, start, rawline))
                    residue_ids.extend(range(start, end+1))
        return residue_ids

    def atom_id_by_residue_id_and_atom_name(self, residue_id, atom_name):
        """For a given residue ID and atom name, return the matching atom ID.
        Return `None` if no match. Exit the program in case of unexpected
        output.
        """
        self._maskstring = ":%s@%s" % (residue_id, atom_name)
        ambmask_stdout = self._run()
        if not ambmask_stdout:
            log.debug("%s did not write stdout." % self)
            return None
        if len(ambmask_stdout.splitlines()) > 1:
            sys.exit("Expected %s to be one line only." % out)
        line = ambmask_stdout.strip()
        if not line.startswith("ATOM"):
            sys.exit("Expected %s to start with ATOM." % line)
        tokens = line.split()
        first = int(tokens[1])
        second = int(tokens [2])
        if first != second:
            sys.exit("Unexpected: identified more than one atom. Output: %s" %
                out)
        return first

    def __repr__(self):
        return ("<%s(topology_file_path=%r, coordinate_file_path=%r, "
            "maskstring=%r, printlevel=%r, outformat=%s)>" % (
            self.__class__.__name__,
            self._topology_file_path,
            self._coordinate_file_path,
            self._maskstring,
            self._printlevel,
            self._outformat))


def restraint_for_residue_id_and_atom_group(resid, atom_group, restraint):
    log.info("Extracting atom IDs for torsional group in residue '%s'" % resid)
    log.info("Atom group: %s" % atom_group)
    log.info("Restraint: %s" % restraint)
    atom_ids = []
    for atom_name in atom_group:
        atom_id = AMBMASK.atom_id_by_residue_id_and_atom_name(resid, atom_name)
        if atom_id is None:
            sys.exit("No atom match for atom %s in group. Exit." % atom_name)
        atom_ids.append(atom_id)
    correspondance = zip(atom_group, atom_ids)
    log.info("Correspondance: %s" % correspondance)
    atom_number_string = ", ".join(str(a) for a in atom_ids)
    iats = "iat=%s" % atom_number_string
    restraint_string="""
&rst
    %s,
    %s,
&end
""" % (iats, restraint)
    return restraint_string


def main():
    global AMBMASK
    try:
        json_input_file = sys.argv[1]
        amber_topology_file = sys.argv[2]
        amber_coordinate_file = sys.argv[3]
    except IndexError:
        sys.stderr.write("Argument error.\n")
        sys.exit("%s\n" % USAGE)
    AMBMASK = Ambmask(amber_topology_file, amber_coordinate_file)

    log.info("Reading '%s'" % json_input_file)
    with codecs.open(json_input_file, encoding="utf-8") as f:
        torsion_configs = json.loads(f.read())
    log.debug("JSON data structure: \n %r" % torsion_configs)

    for config in torsion_configs:
        for residue_name in config["residue_names"]:
            log.info("Extracting residue IDs for residue '%s'" % residue_name)
            resids = AMBMASK.residue_ids_by_name(residue_name)
            if resids is None:
                log.info("No such residue in the system. Skip.")
                continue
            log.info("Residue IDs: %s" % ", ".join("%s" % r for r in resids))
            for resid in resids:
                for (atom_group, restraint) in config["torsion_groups"]:
                    log.info("Building restraint string for group %s" %
                        atom_group)
                    restraint_string = restraint_for_residue_id_and_atom_group(
                        resid,
                        atom_group,
                        restraint
                        )
                    print restraint_string


AMBMASK = None
USAGE="""
Usage:
%s json_config_file amber_topology_file amber_coordinate_file
""" % os.path.basename(sys.argv[0])




# The toplevel config values contain dihedral names. A dihedral name must
# consist of the names of the residues involved in a dihedral. If a dihedral
# involves multiple residues, the corresponding residue names must be linked
# with a dash and written in the right order.
example_config="""
04V-4ZB:
    psi:
        resnames: [04V, 04V, 4ZB, 4ZB]
        atoms: [H1, C1, O4, C4]
    phi:
        resnames: [04V, 4ZB, 4ZB, 4ZB]
        atoms: [H1, C1, O4, C4]
4ZB-34V:
    psi:
        resnames: [4ZB, 4ZB, 34V, 34V]
        atoms: [H1, C1, O3, C3]
    phi:
        resnames: [4ZB, 34V, 34V, 34V]
        atoms: [C1, O3, C3, H3]
34V-4ZB:
    psi:
        resnames: [34V, 34V, 4ZB, 4ZB]
        atoms: [H1, C1, O4, C4]
    phi:
        resnames: [34V, 4ZB, 4ZB, 4ZB]
        atoms: [C1, O4, C4, H4]
"""

ambmask = None

def main2():
    global ambmask
    parser = argparse.ArgumentParser(description='Dickes Tool.')
    parser.add_argument('topologyfile', action="store")
    parser.add_argument('coordinatefile', action="store")
    parser.add_argument('trajectoryfile', action="store")
    parser.add_argument('configfile', action="store")
    parser.add_argument('-i', '--inverse-search',
        action='store_true', default=False)
    options = parser.parse_args()
    print options
    ambmask = Ambmask(
        topology_file_path=options.topologyfile,
        coordinate_file_path=options.coordinatefile)



    with codecs.open(options.configfile, encoding="utf-8") as f:
        c = f.read()
        log.info("Read configuration file: \n%s" % c)
        config = yaml.load(c)
    validate_config(config)
    resname_resids_mapping = get_resids_for_resnames(config)
    log.info("Residue IDs as identified by ambmask: %s" %
        resname_resids_mapping)
    dihedrals = identify_dihedrals(
        config,
        resname_resids_mapping,
        options.inverse_search)
    dihedral_prettystring = "\n".join(str(d) for d in dihedrals)
    log.info("Identified %s dihedrals: \n%s" % (
        len(dihedrals), dihedral_prettystring))


class Atom(object):
    def __init__(self, aid, aname, resid, resname):
        self.aid = aid
        self.name = aname
        self.resid = resid
        self.resname = resname

    def __str__(self):
        return "%s_%s" % (str(self.aid).zfill(5), self.name)

    def __repr__(self):
        return "Atom(aid=%s, aname=%s, resid=%s, resname=%s)" % (
            self.aid, self.name, self.resid, self.resname)


class Dihedral(object):
    def __init__(self, name, atoms):
        self.name = name
        assert len(atoms) == 4
        for a in atoms:
            assert isinstance(a, Atom)
        self.atoms = atoms

    def __str__(self):
        atomstring = "--".join(str(a) for a in self.atoms)
        return "dihedral %s with atoms %s" % (self.name, atomstring)

    def __repr__(self):
        return "Dihedral(name=%s, atoms=%s)" % (
            self.name, self.atoms)


def identify_dihedrals(config, resname_resids_mapping, inverse):
    """For each dihedral name as given in the config, find one or multiple
    residue units, i.e. residues that are in sequence (as of their residue ID)
    and correspond to the residue names as defined in the dihedral name. For
    instance, for a dihedral with name "RUX-RAX", find all occurences of
    neighboring RUX and RAX residues, whereas the order must be maintained
    in search direction.

    For each residue unit corresponding to one dihedral name, build a
    representation of the specific dihedral containing all important details,
    especially atomic IDs that are valid in context of the Amber topology file.
    Based on the residue IDs within one torsion residue unit and the
    corresponding atom names as specified in the config, Ambpdb is used to
    retrieve the atomic IDs of the atoms involved in one specific dihedral.
    An abstract dihedral object is created containing all information about
    this dihedral such as involved residue names, residue numbers, atom names,
    and atom numbers.

    Return a list of dihedral objects.

    This method supposes and requires that the residues as defined in a
    dihedral name must be in sequence, i.e. that the corresponding residue
    IDs (integer numbers) numbers must not have gaps inbetween.

    The search for units by default goes from large to small residue numbers
    (for GAGs, the nonreducing to reducing end direction usually goes from
    large to small, since ROH usually is the first residue).
    The search direction can be inversed via `inverse`.
    """
    # Create a list of potentially involved residues, whereas each element is
    # 2-tuple and contains residue name as well as ID, such as ('34V', 11).
    resname_id_tuples = []
    for resname, ids in resname_resids_mapping.iteritems():
        resname_id_tuples.extend([(resname, idd) for idd in ids])
    # `resname_id_tuples` now is something like
    # [('04V', 13), ('34V', 3), ('34V', 5), ('34V', 7), ('34V', 9), ...]
    # Sort this list by residue ID.
    # If `inverse` is False (default), sorting should happen from large to
    # small, i.e. in reversed direction from Python's point of view.
    resname_id_tuples.sort(key=operator.itemgetter(1), reverse=not inverse)
    log.info("Sorted residue list in dihedral unit search direction:\n%s" %
        resname_id_tuples)
    log.info(("Note: specify the --inverse-search option in order to reverse "
        "the search direction."))

    # Find dihedral residue units in this list by moving a selective window
    # over it. The witdth of the window corresponds to the number of residues
    # involved in a dihedral type. In case a match is found, validate it and
    # pull missing details (atom IDs) from the topology file via Ambpdb. Create
    # a `Dihedral` object and store it.
    dihedrals = []
    for dihedname in config.keys():
        involved_residue_names = tuple(dihedname.split('-'))
        log.info("Looking for matches for dihedral %s." % dihedname)
        all_units_of_right_length = window(resname_id_tuples,
            len(involved_residue_names))
        for potential_unit in all_units_of_right_length:
            unit_resnames, unit_resids = zip(*potential_unit)
            if unit_resnames == involved_residue_names:
                log.info("Unit matches: %s." % (potential_unit, ))
                # Validate that the residue IDs are in sequence:
                unit_resids_ascending = sorted(unit_resids)
                log.debug("Residue IDs involved in unit (ascending): %s." %
                    unit_resids_ascending)
                # http://stackoverflow.com/q/18131741/145400
                assert unit_resids_ascending == range(unit_resids_ascending[0],
                    unit_resids_ascending[-1]+1), (("Identified dihedral "
                        "residue unit with non-coherent sequence IDs: %s.") %
                        unit_resids_ascending)
                # Create Dihedral object.
                # Each dihedral unit might involve multiple angles, two in
                # this case, phi and phi:
                #
                # 04V-4ZB:
                #     psi:
                #         resnames: [04V, 04V, 4ZB, 4ZB]
                #         atoms: [H1, C1, O4, C4]
                #     phi:
                #         resnames: [04V, 4ZB, 4ZB, 4ZB]
                #         atoms: [H1, C1, O4, C4]

                # Create dictionary from potential unit, enabling easy lookup
                # of resid for given resname in current unit.
                unit = dict(potential_unit)
                dihed_config = config[dihedname]
                for angle_name, details in dihed_config.iteritems():
                    # Extract / generate useful lists/tuples.
                    atomnames = details['atoms']
                    resnames = details['resnames']
                    resids = [unit[n] for n in resnames]
                    involved_resids = [unit[n] for n in involved_residue_names]
                    residue_atom_name_groups = zip(resids, resnames, atomnames)
                    # Create descriptive label for current dihedral, containing
                    # the two residues (incl. id and name) as well as the name
                    # of the angle as specified in the config.
                    involved_residues_identifier = "-".join(
                        "%s_%s" % (str(i).zfill(3),n) for i,n in zip(
                        involved_resids, involved_residue_names))
                    dihedral_identifier = "%s-%s" % (
                        involved_residues_identifier, angle_name)
                    log.info("Creating dihedral %s involving atoms %s." % (
                        dihedral_identifier, atomnames))
                    log.info("Atom name - residue matching: %s." %
                        residue_atom_name_groups)
                    #log.debug("Current dihedrals's ")
                    atoms = []
                    log.info("Looking up atom IDs via ambmask...")
                    for resid, resname, aname in residue_atom_name_groups:
                        aid = ambmask.atom_id_by_residue_id_and_atom_name(
                            resid, aname)
                        a = Atom (aid, aname, resid, resname)
                        atoms.append(a)
                        log.debug("Created atom %r." % a)
                    d = Dihedral(
                        name=dihedral_identifier,
                        atoms=atoms)
                    log.info("Created dihedral object %r." % d)
                    dihedrals.append(d)
                    #break
                #break
            #break
    return dihedrals




def get_resids_for_resnames(config):
    """Return a dictionary containing the residue IDs for every residue
    ocurring in any dihedral.
    """
    dihedral_names = config.keys()
    resnames_in_dihedral_names = set(itertools.chain.from_iterable(
        s.split('-') for s in dihedral_names))
    log.info("Residues in dihedral names: %s" % resnames_in_dihedral_names)
    resname_resids_mapping = {}
    for name in resnames_in_dihedral_names:
        resname_resids_mapping[name] = ambmask.residue_ids_by_name(name)
    return resname_resids_mapping


def validate_config(config):
    for dihedname in config:
        resnames = []
        for angle in  config[dihedname]:
            assert 'atoms' in config[dihedname][angle]
            assert 'resnames' in config[dihedname][angle]
            assert len(config[dihedname][angle]['resnames']) == 4
            assert len(config[dihedname][angle]['atoms']) == 4
            resnames.extend(config[dihedname][angle]['resnames'])
        resnames_in_details = set(resnames)
        resnames_in_dihedral_name = set(dihedname.split('-'))
        if not resnames_in_details == resnames_in_dihedral_name:
            sys.exit(("Dihedral name %s contains other residues than specified"
                " in dihedral name: %s" % (dihedname, resnames_in_details)))
        log.debug(
            "Residue names in details for dihedral %s match." % dihedname)


def window(iterable, n):
    # Based on http://stackoverflow.com/a/7636587/145400
    els = itertools.tee(iterable, n)
    for i, el in enumerate(els):
        for _ in xrange(i):
            next(el, None)
    return itertools.izip(*els)


if __name__ == "__main__":
    main2()


