import os
import subprocess
import gzip
import shutil

class DatabaseDownloader:
    """Acts as a handler for HMM and sequence databases to search FASTAs against.
    """
    def __init__(self, output_dir, remove_zips=False):
        """
        Args:
            output_dir (str): Path to the directory you want to put databases files in
            remove_zips (bool): Option to delete compressed files and folders after unpacking to correct places
        """
        self.output_dir = output_dir
        self.remove_zips = remove_zips
        self.dbs_to_dl = {"uniref50.fasta.gz": "ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz",
                    "profiles.tar.gz": "ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz",
                    "ko_list.gz": "ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz",
                    "Pfam-A.hmm.gz": "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz",
                    "Pfam-A.hmm.dat.gz": "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.dat.gz",
                    "NCBIfam-AMRFinder.HMM.tar.gz": "https://ftp.ncbi.nlm.nih.gov/hmm/NCBIfam-AMRFinder/latest/NCBIfam-AMRFinder.HMM.tar.gz",
                    "vog.hmm.tar.gz": "http://fileshare.csb.univie.ac.at/vog/latest/vog.hmm.tar.gz",
                    "bacteria_odb10.2019-06-26.tar.gz": "https://busco-data.ezlab.org/v4/data/lineages/bacteria_odb10.2019-06-26.tar.gz",
                    "archaea_odb10.2019-01-04.tar.gz": "https://busco-data.ezlab.org/v4/data/lineages/archaea_odb10.2019-01-04.tar.gz",
                    "eukaryota_odb10.2019-11-20.tar.gz": "https://busco-data.ezlab.org/v4/data/lineages/eukaryota_odb10.2019-11-20.tar.gz"}

    def check_if_downloaded(self):
        """Checks if all databases in self.dbs_to_dl have already been downloaded

        Returns:
            not_downloaded (dict): Databases that have not been downloaded yet, along with their respective links
        """
        not_downloaded = {}
        for db_name, db_url in self.dbs_to_dl.items():
            if not os.path.exists(os.path.join(self.output_dir, db_name)):
                not_downloaded[db_name] = db_url
        return not_downloaded

    def download_file(self, url, output_filename):
        download_cmd = ['wget', '-O', os.path.join(self.output_dir, output_filename), url]
        subprocess.run(download_cmd)

    def retrieve_databases(self):
        """Checks whether databases have been downloaded, then downloads ones that haven't
        """
        files_to_dl = self.check_if_downloaded()
        for db_name, db_url in files_to_dl.items():
            print('Downloading ' + db_name)
            self.download_file(db_url, db_name)

    def remove_file(self, file_loc):
        if self.remove_zips:
            os.remove(file_loc)

    def convert_uniref_to_mmseqs(self):
        """Converts UniRef50 database to MMSeqs2 db
        """
        print('Converting UniRef50 to MMSeqs2 db')
        uniref_mmseqs_loc = os.path.join(self.output_dir, "uniref_mmseqs.db")
        if not os.path.exists(uniref_mmseqs_loc):
            uniref_full_loc = os.path.join(self.output_dir, "uniref50.fasta.gz")
            createdb_cmd = ['mmseqs', 'createdb', uniref_full_loc, uniref_mmseqs_loc]
            subprocess.run(createdb_cmd)
            self.remove_file(uniref_full_loc)

    def unpack_kofam_db(self):
        """Unzips the KOfam HMM tarball and writes individual HMMs to master HMM file
        """
        kofam_concatenated_hmms_loc = os.path.join(self.output_dir, 'kofam.hmm')
        if not os.path.exists(kofam_concatenated_hmms_loc):
            kofam_loc = os.path.join(self.output_dir, "profiles.tar.gz")
            kofam_profile_dir = os.path.join(self.output_dir, "kofam_profiles/")
            if not os.path.exists(kofam_profile_dir):
                os.mkdir(kofam_profile_dir)
                unpack_cmd = ['tar', '-xvzf', kofam_loc, '-C', kofam_profile_dir]
                subprocess.run(unpack_cmd)
            print('Unpacking KOfam HMMs')
            kofam_profiles_subdir = os.path.join(kofam_profile_dir, "profiles")
            with open(kofam_concatenated_hmms_loc, 'w') as kofam:
                for filename in os.listdir(kofam_profiles_subdir):
                    if filename.endswith('.hmm'):
                        with open(os.path.join(kofam_profiles_subdir, filename), 'r') as hmm:
                            kofam.write(hmm.read())
            self.remove_file(kofam_loc)
            self.remove_file(kofam_profile_dir)

    def unpack_kofam_list(self):
        """Unzips the KOfam HMM description list
        This file has the specific e-value and bit score cutoffs for each HMM in the KOfam list.
        """
        print('Unpacking KOfam gene list')
        kofam_list_output_loc = os.path.join(self.output_dir, "ko_list.tsv")
        if not os.path.exists(kofam_list_output_loc):
            kofam_list_loc = os.path.join(self.output_dir, "ko_list.gz")
            with gzip.open(kofam_list_loc, 'rb') as list_in:
                with open(kofam_list_output_loc, 'wb') as list_out:
                    shutil.copyfileobj(list_in, list_out)
            self.remove_file(kofam_list_loc)

    def unpack_pfam_db(self):
        """Unzips the Pfam-A HMM tarball
        """
        pfam_concatenated_hmms_loc = os.path.join(self.output_dir, "pfam.hmm")
        if not os.path.exists(pfam_concatenated_hmms_loc):
            pfam_loc = os.path.join(self.output_dir, "Pfam-A.hmm.gz")
            with gzip.open(pfam_loc, 'rb') as pfam_in:
                with open(pfam_concatenated_hmms_loc, 'wb') as pfam_out:
                    shutil.copyfileobj(pfam_in, pfam_out )
            self.remove_file(pfam_loc)

    def unpack_amrfinder_db(self):
        """Unzips the AMRfinder HMM tarball and writes individual HMMs to master HMM file
        """
        print('Unpacking AMRfinder db')
        amrfinder_concatenated_hmms_loc = os.path.join(self.output_dir, "amrfinder.hmm")
        if not os.path.exists(amrfinder_concatenated_hmms_loc):
            amrfinder_loc = os.path.join(self.output_dir, "NCBIfam-AMRFinder.HMM.tar.gz")
            amrfinder_profile_dir = os.path.join(self.output_dir, "amrfinder_profiles/")
            if not os.path.exists(amrfinder_profile_dir):
                os.mkdir(amrfinder_profile_dir)
                unpack_cmd = ['tar', '-xvzf', amrfinder_loc, '-C', amrfinder_profile_dir]
                subprocess.run(unpack_cmd)
            amrfinder_profiles_subdir = os.path.join(amrfinder_profile_dir, 'HMM')
            with open(amrfinder_concatenated_hmms_loc, 'w') as amrfinder:
                for filename in os.listdir(amrfinder_profiles_subdir):
                    if filename.endswith('.HMM'):
                        with open(os.path.join(amrfinder_profiles_subdir, filename), 'r') as hmm:
                            amrfinder.write(hmm.read())
            self.remove_file(amrfinder_loc)
            self.remove_file(amrfinder_profile_dir)

    def unpack_vog_db(self):
        """Unzips the VOGdb HMM tarball and writes individual HMMs to master HMM file
        """
        print('Unpacking VOGdb')
        vog_concatenated_hmms_loc = os.path.join(self.output_dir, "vog.hmm")
        if not os.path.exists(vog_concatenated_hmms_loc):
            vog_loc = os.path.join(self.output_dir, "vog.hmm.tar.gz")
            vog_profile_dir = os.path.join(self.output_dir, "vog_profiles/")
            if not os.path.exists(vog_profile_dir):
                os.mkdir(vog_profile_dir)
                unpack_cmd = ['tar', '-xvzf', vog_loc, '-C', vog_profile_dir]
                subprocess.run(unpack_cmd)
            with open(vog_concatenated_hmms_loc, 'w') as vog:
                for filename in os.listdir(vog_profile_dir):
                    if filename.endswith('.hmm'):
                        with open(os.path.join(vog_profile_dir, filename), 'r') as hmm:
                            vog.write(hmm.read())
            self.remove_file(vog_profile_dir)

    def unpack_busco_dbs(self):
        '''Unzips Bacteria, Archaea, and Eukaryota BUSCO dbs and places them in a folder'''
        print('Unpacking BUSCO databases')
        busco_dir = os.path.join(self.output_dir, "busco_databases")
        bac_loc = os.path.join(self.output_dir, "bacteria_odb10.2019-06-26.tar.gz")
        arc_loc = os.path.join(self.output_dir, "archaea_odb10.2019-01-04.tar.gz")
        euk_loc = os.path.join(self.output_dir, "eukaryota_odb10.2019-11-20.tar.gz")
        if not os.path.exists(busco_dir):
            os.mkdir(busco_dir)
        bac_dir = os.path.join(busco_dir, "bacteria")
        if not os.path.exists(bac_dir):
            os.mkdir(bac_dir)
            unpack_cmd = ['tar', '-xvzf', bac_loc, '-C', bac_dir]
            subprocess.run(unpack_cmd)
        arc_dir = os.path.join(busco_dir, "archaea")
        if not os.path.exists(arc_dir):
            os.mkdir(arc_dir)
            unpack_cmd = ['tar', '-xvzf', arc_loc, '-C', arc_dir]
            subprocess.run(unpack_cmd)
        euk_dir = os.path.join(busco_dir, "eukaryota")
        if not os.path.exists(euk_dir):
            os.mkdir(euk_dir)
            unpack_cmd = ['tar', '-xvzf', euk_loc, '-C', euk_dir]
            subprocess.run(unpack_cmd)
        self.remove_file(bac_loc)
        self.remove_file(arc_loc)
        self.remove_file(euk_loc)

    def download_and_unpack_databases(self):
        """Workhorse function that performs retrieval and unpacking steps for databases
        """
        self.retrieve_databases()
        self.unpack_kofam_db()
        self.unpack_kofam_list()
        self.unpack_pfam_db()
        self.unpack_amrfinder_db()
        self.unpack_vog_db()
        self.convert_uniref_to_mmseqs()