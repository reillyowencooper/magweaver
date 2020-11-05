import argparse, os, sys, logging
import src.mag_decontaminate as decon
import src.database_handling as db

__author__ = "Reilly Cooper"
__maintainer__ = "Reilly Cooper"
__email__ = "reilly.owen.cooper@gmail.com"
__version__ = "0.0.1"

def help_printer():
    print('')
    print('    ~~~~ CAMAG v' + __version__ + ' ~~~~    ')
    print('')
    print('''\
        
        Choose a workflow or a module (currently only selection and decontamination is working).
        To get help for a workflow, add the -h flag; for example, camag decontaminate -h
        
        Workflows:
            cmag: (END GOAL) Evaluate genome bins for completeness, then refine those bins for attempted CMAGs. 
            annotate: (END GOAL): Annotate genomes based on their taxonomic identity
        
        Modules:  
            decontaminate: Remove contigs with unlikely metrics from a MAG
            quality: Assess genome quality (completeness/contamination) using CheckM
        
        ''')
    print('')
    print('    ~~~~ Good luck! ~~~~        ')
    print('')
    
def parse_args(args):
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='ops')
    parent_parser = argparse.ArgumentParser()
    
    # Decontaminate parent
    decontaminate_parent = argparse.ArgumentParser(add_help = False)
    decon_flags = decontaminate_parent.add_argument_group("Decontamination options")
    # Will add arguments corresponding to deviation from mean cov/gc/tetra, probably add option to change taxonomic rank level
    decon_flags.add_argument('-s', '--suspicion', 
                             help = 'Contig suspicion level to flag corresponding to number of metrics outside of mean. Min 1, max 4', 
                             default = 3, type = int)
    
    # Quality parent
    quality_parent = argparse.ArgumentParser(add_help = False)
    quality_flags = quality_parent.add_argument_group("Quality options")
    quality_flags.add_argument('-comp', '--completion', 
                               help = 'Minimum completion to consider genome highquality', 
                               default = 95.0, type = float)
    quality_flags.add_argument('-contam', '--contamination',
                               help = 'Maximum contamination (redundancy) to consider genome high quality',
                               default = 5.0, type = float)
    quality_flags.add_argument('-contig', '--num_contigs',
                               help = 'Maximum number of contigs allowed in genome to consider high quality',
                               default = 10, type = int)
    
    # Decontaminate module parser
    decontaminate_parser = subparsers.add_parser("decontaminate", 
                                                 parents = [parent_parser, decontaminate_parent],
                                                 add_help = False)
    
    decon_io = decontaminate_parser.add_argument_group("Input/Output")
    decon_io.add_argument('-i', '--input_mag', help = 'Path to input MAG', type = str)
    decon_io.add_argument('-fwd', '--forward_reads', help = 'Path to forward reads used for assembly', type = str)
    decon_io.add_argument('-rev', '--reverse_reads', help = 'Path to reverse reads used for assembly', type = str)
    decon_io.add_argument('-o', '--output_dir', help = 'Directory to output resulting files', type = str)
    
    # Quality module parser
    quality_parser = subparsers.add_parser("quality",
                                           parents = [parent_parser, quality_parent],
                                           add_help = False)
    
    quality_io = quality_parser.add_argument_group("Input/Output")
    quality_io.add_argument('-i', '--input_dir', help = 'Directory with MAGs to analyze', type = str)
    quality_io.add_argument('-o', '--output_dir', help = 'Directory to output resulting files', type = str)
    
    # Parse user input
    if len(args) < 1 or args[0] == '-h':
        help_printer()
        sys.exit()
    else:
        return parser.parse_args(args)
    
class Shuttler(object):
    '''Take argparse inputs and enacts correct commands'''
    def __init__(self, swissprot_db, scg_db):
        self.logger = logging.getLogger()
        self.swissprot_db = swissprot_db
        self.scg_db = scg_db
        
    def create_logger(self, **kwargs):
        self.logger = logging.getLogger()
        #
        #logging.basicConfig(format = '%(asctime)s %(levelname)-8s %(message)s',
        #                    level = logging.DEBUG,
        #                    filename = os.path.join(kwargs['output_dir'], "camag.log"))
        #logging.debug('CAMAG started running with cmd: {0}'.format(' ').join(**kwargs))
        
        terminal = logging.StreamHandler()
        terminal.setLevel(logging.INFO)
        terminal_format = logging.Formatter('%(asctime)s %(message)s')
        terminal.setFormatter(terminal_format)
        logging.getLogger('').addHandler(terminal)
                  
    def quality_ops(self, **kwargs):
        logging.debug('Analyzing MAG bin quality')
        ### Do something ###
        logging.debug('Finished analyzing MAG quality, high quality bins deposited in ' + kwargs['output_dir'])
        
    def decontam_ops(self, **kwargs):
        logging.debug('Decontaminating ' + kwargs['input_mag'])
        decon.decontaminate_wrapper(self.swissprot_db, self.scg_db, **kwargs)
        logging.debug('Completed decontamination of ' + kwargs['input_mag'])
        
    def route_args(self, args):
        self.create_logger()
        if args.ops == 'quality':
            self.quality_ops(**vars(args))
        elif args.ops == 'decontaminate':
            self.decontam_ops(**vars(args))
            
            