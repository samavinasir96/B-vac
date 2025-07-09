import pandas as pd
from Bio import SeqIO
from tkinter import *
from tkinter import messagebox
from tkinter import ttk
import os
from time import sleep
from io import StringIO
from Bio import Align
from ttkwidgets.autocomplete import AutocompleteCombobox
from tkinter import filedialog
import threading
import swifter
import time

start_time = time.time()

window = Tk()
window.title("B-Vac")
window.size()
window.geometry("1100x660")
window.maxsize(height=660, width=1100)

Icon_path = PhotoImage(file="./Images/Logo/Email.png")
window.iconphoto(False,Icon_path)
HELP_ICON = PhotoImage(file="./Images/Logo/Help.png")
SAVE_BUTTON = PhotoImage(file="./Images/Logo/Save_Button.png")
BROWSE_BUTTON = PhotoImage(file="./Images/Logo/Browse_Button.png")
SUBMIT_BUTTON = PhotoImage(file="./Images/Logo/SUBMIT.png")
IBG_LAB_LOGO = PhotoImage(file="./Images/New IBG Logo.png")
MGBio_LOGO = PhotoImage(file="./Images/MGBio.png")

aligner = Align.PairwiseAligner()

First_Color = Second_Color = "#EAEAEA"
Buttons_Colors = Label_Color = "#003380"
Button_Text = "#F5F5F5"
font_family = "Bahnschrift"
width_frames = 1150

Submission_method_frame = Frame(window, width=width_frames, height=250, bg=First_Color)
Submission_method_frame.grid(row=0,sticky="w")

IBG_LOGO = Label(Submission_method_frame, image=IBG_LAB_LOGO, bg=Second_Color)
IBG_LOGO.place(relx=0.855, rely=0.09)

MGBio_LOGO_Label = Label(Submission_method_frame, image=MGBio_LOGO, bg=Second_Color)
MGBio_LOGO_Label.place(relx=0.79, rely=0.1)

Job_Status = Label(Submission_method_frame, text="Job Status: Not Running",bg=First_Color, font=(font_family,(14)), fg=Buttons_Colors)
Job_Status.place(relx=0.39, rely=0.867)

Logo_Label_Image = PhotoImage(file="./Images/Logo/Logo_Final.png")
Logo_Label = Label(Submission_method_frame, image=Logo_Label_Image, bg=First_Color)
Logo_Label.place(relx=0.022, rely=0.07)

Description_Label = Label(Submission_method_frame, text="A Robust Software Package for Bacterial Vaccine Design", bg=First_Color, font=(font_family,(16)), fg=Buttons_Colors)
Description_Label.place(relx=0.25, rely=0.16)

Submission_label = Label(Submission_method_frame, text="Choose a Submission Method", bg=First_Color, font=(font_family,(10)), fg=Buttons_Colors)
Submission_label.place(relx=0.022, rely=0.3)


def all_equal(items):
    """Returns True iff all items are equal."""
    first = items[0]
    return all(x == first for x in items)

def compute_matched(aligned_sequences):
    match_count = sum(1 for chars in zip(*aligned_sequences) if all_equal(chars))
    total = len(aligned_sequences[0])
    mismatch_count = total - match_count  # Obviously.
    percentage = int(match_count / total * 100)
    return percentage

def read_submitted_fasta_file(file_name):
    title_list = []
    fasta_list = []
    description_list = []
    for seq_record in SeqIO.parse(file_name, "fasta"):
        title_list.append(seq_record.id)
        fasta_list.append(seq_record.seq)
        description_list.append(seq_record.description)
    return title_list,fasta_list, description_list

def read_submitted_fasta_sequence(sequence):
    title_list = []
    fasta_list = []
    fasta_io = StringIO(sequence)
    records = SeqIO.parse(fasta_io, "fasta")
    for record in records:
        title_list.append(record.id)
        fasta_list.append(str(record.seq))
    return title_list,fasta_list

string_fasta = """>WP_003355768.1
MRLRKKWWARPEIEASDKFAEEPKELRGKWNKEFNNNNDIHLELGCGRGGFISQLVEKNKDINYVGIDLKDEVIVYAIRK
VKEKEEEVKREFKNIKFVTMNIMGIAEVFDKNEISKIYINFCNPWPKERHNKRRLTHTKLLTEYKKFLKPNTEIWFKTDD
KELFEDSQEYFKESGFNIEYITYDLHNSDFKENIKTEYETKFETMGMKIMFLKARLL\n"""

fasta_sequence = read_submitted_fasta_sequence(string_fasta)

Fasta_Entry_Var = StringVar()
Fasta_Entry = Text(Submission_method_frame, height=5.5, width=105,padx=8, pady=8)
Fasta_Entry.place(rely=0.4, relx=0.025)

readed_file_data = []

filename_var = StringVar()

def FASTA_BROWSE():

    try:
        filename = list(filedialog.askopenfilenames(parent=window, initialdir="/Desktop", title="Choose Fasta File", filetypes=[("fasta file(*.FAA)", "*.FAA"), ("fasta file(*.FASTA)", "*.FASTA")]))

        Total_sequences_length = []
        for files in filename:
            Fasta_Entry.delete('1.0', END)
            fasta_sequences_title, fasta_sequences_file, fasta_sequences_des = read_submitted_fasta_file(files)
            Total_sequences_length.append(len(fasta_sequences_file))

        FASTA_BROWSE_LABEL.config(text=f"{len(Total_sequences_length)} File(s) >{sum(Total_sequences_length)}")
        filename_var.set(str(filename))

    except AttributeError:
        messagebox.showwarning("Warning", "No File Selected")


Fasta_Browse = Button(Submission_method_frame, image=BROWSE_BUTTON, command=FASTA_BROWSE, bd=0)
Fasta_Browse.place(relx=0.79, rely=0.4)

FASTA_BROWSE_LABEL = Label(Submission_method_frame, text="Select File", height=2, width=20, bd=0, background=First_Color, font = (font_family,12), fg=Buttons_Colors)
FASTA_BROWSE_LABEL.place(relx=0.78, rely=0.549)


def SET_DEFAULT_VALUES():
    Localization_identity.set(60)
    Reliability.set(50)
    Non_Homologs_identity.set(35)
    Non_Homologs_perc.set(50)
    VFDB_identity.set(70)
    Epitope_identity_CD4.set(40)
    Epitope_identity_CD8.set(50)
    Epitope_window_CD4.set(20)
    Epitope_window_CD8.set(9)
    Epitope_window_B_cell.set(15)
    Epitope_identity_B_cell.set(50)

def FILL_EXAMPLE():
    Fasta_Entry.insert(END, string_fasta)

def MUST_EVALUATE():
    if Must_Evalutate_var.get() == "True":
        messagebox.showinfo("Info", "Must Evaluate Selected")
        Localization.select()
        Non_Homologs_check_btn.select()
        VFDB_check_btn.select()
        Epitope_check_btn.select()
        Localization_combo_box_1_var.set("All")

    else:
        Localization.deselect()
        Non_Homologs_check_btn.deselect()
        VFDB_check_btn.deselect()
        Epitope_check_btn.deselect()


SAVE_FILE_Var = StringVar()

def SAVE_FILE():
    folder_selected = filedialog.askdirectory()
    SAVE_FILE_Var.set(folder_selected)


def Update_Progress(value):
    my_progress["value"] = value
    window.update()

def FASTA_SUBMIT_BTN():

    data_list = []

    input = Fasta_Entry.get("1.0",END).upper()
    if len(input) == 1 and filename_var.get() != "":

        all_names = filename_var.get().strip('][').replace("'", "").split(", ")

        all_titles = []
        all_fastas = []

        for name in all_names:
            titles, fastas, desc = read_submitted_fasta_file(name)
            for title, fasta in zip(titles, fastas):
                all_titles.append(title)
                all_fastas.append(fasta)

        data_list.append({
            "TITLE": all_titles,
            "FASTA": all_fastas
        })

        FASTA_BROWSE_LABEL.config(text="Select File")

    elif len(input) != 0 and filename_var.get() != "":
        messagebox.showwarning("Warning", "Choose One Submission Method")
        filename_var.set("")
        FASTA_BROWSE_LABEL.config(text="Select File")

    else:
        Amino_acid_list = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y",
                           "V"]
        fasta_io = StringIO(input)
        records = SeqIO.parse(fasta_io, "fasta")
        sure_var = ""
        for record in records:
            for single_amino_input in list(record.seq):
                if single_amino_input in Amino_acid_list:
                    sure_var = True
                else:
                    sure_var = False
                    break
            if sure_var == True:
                title, fasta = read_submitted_fasta_sequence(input)
                data_list.append({
                    "TITLE" : title,
                    "FASTA" : fasta
                })
                break

            else:
                messagebox.showerror("Error", "Mistake in Fasta Sequence")
                Fasta_Entry.delete('1.0', END)

    if len(data_list) == 0 and filename_var.get() == "":
        messagebox.showerror("Warning", "Choose a Submission Method")
    else:

        if Must_Evalutate_var.get() != "True":

            Job_Status.config(text="Job Status: Running", fg="Red")

            if SAVE_FILE_Var.get() != "":

                Progress_bar_list = []

                if Localization_var.get() != "True" and Non_Homologs_var.get() != "True" and VFDB_var.get() != "True" and Epitope_var.get() != "True":
                    messagebox.showerror("Error", "Choose A Method")

                if Localization_var.get() == "True":
                    Progress_bar_list.append(25)

                if Non_Homologs_var.get() == "True":
                    Progress_bar_list.append(25)

                if VFDB_var.get() == "True":
                    Progress_bar_list.append(25)

                if Epitope_var.get() == "True":
                    Progress_bar_list.append(25)


                Bar_Value = []

                try:
                    Bar = 100/len(Progress_bar_list)
                    Bar_Value.append(Bar)
                except ZeroDivisionError:
                    Bar_Value.append(100)

                Update_Progress(Bar_Value[0])

                if Localization_var.get() == "True":

                    LOCATION_LIST = []
                    BIT_SCORE_LIST = []
                    PROTEIN_ID_LIST = []
                    FASTA_LIST = []

                    TITLE, FASTA, DESCRIPTION = read_submitted_fasta_file(
                        f"./Data/FINAL_FILES/{FILE_NAME_VAR.get()}")

                    for TITLE_READED, FASTA_READED, DESCRIPTION_READED in zip(TITLE, FASTA, DESCRIPTION):
                        info = " ".join(DESCRIPTION_READED.split(" ")[1:-1])
                        score = DESCRIPTION_READED.split(" ")[-1]

                        LOCATION_LIST.append(info)
                        BIT_SCORE_LIST.append(int(score))
                        PROTEIN_ID_LIST.append(TITLE_READED)
                        FASTA_LIST.append(FASTA_READED)

                    LOCATION_DATA = pd.DataFrame({
                        "PROTEIN ID": PROTEIN_ID_LIST,
                        "BIT SCORE": BIT_SCORE_LIST,
                        "LOCATION": LOCATION_LIST,
                        "FASTA": FASTA_LIST
                    }).reset_index(drop=True)

                    if Reliability_var.get() != 0:
                        mean_filtered = LOCATION_DATA[LOCATION_DATA["BIT SCORE"] > Reliability_var.get()]
                        remove_list = ["secreted", "outer membrane", "fimbrium"]
                        NEW_DATA_FRAME = mean_filtered[mean_filtered["LOCATION"].isin(remove_list)].reset_index(
                            drop=True)
                        REALIABILITY_VAR.set(Reliability_var.get())
                        DATA_FINAL_LOC = pd.DataFrame({})
                        LIST_LOC = []

                        for single_Localization in data_list:
                            for single_Localization_fasta, single_Localization_id in zip(single_Localization["FASTA"], single_Localization["TITLE"]):
                                DATA_FINAL_LOC["LOC"] = NEW_DATA_FRAME["LOCATION"]
                                DATA_FINAL_LOC["FASTA_UNIPROT"] = NEW_DATA_FRAME["FASTA"]
                                DATA_FINAL_LOC['FASTA_SUBMITTED'] = single_Localization_fasta
                                DATA_FINAL_LOC['FASTA ID'] = single_Localization_id
                                DATA_FINAL_LOC["LEN_UNIPROT"] = DATA_FINAL_LOC["FASTA_UNIPROT"].apply(len)
                                DATA_FINAL_LOC["LEN_SUBMITTED"] = DATA_FINAL_LOC['FASTA_SUBMITTED'].apply(len)
                                FILTERED_LENGTH = DATA_FINAL_LOC[DATA_FINAL_LOC["LEN_SUBMITTED"] <= DATA_FINAL_LOC["LEN_UNIPROT"]].reset_index(drop=True)
                                try:
                                    FILTERED_LENGTH['PERCENTAGE'] = FILTERED_LENGTH.apply(lambda x: compute_matched([x["FASTA_UNIPROT"], x['FASTA_SUBMITTED']]), axis=1)
                                    FINAL = FILTERED_LENGTH[FILTERED_LENGTH['PERCENTAGE'] == FILTERED_LENGTH['PERCENTAGE'].max()]
                                    FINAL_SEP = FINAL[FINAL["LEN_UNIPROT"] == FINAL[["LEN_UNIPROT", "LEN_SUBMITTED"]]["LEN_UNIPROT"].min()].reset_index(drop=True)
                                    LIST_LOC.append(FINAL_SEP)
                                except ValueError:
                                    pass
                        FINAL_LOC = pd.concat(LIST_LOC)

                        if Localization_identity_var.get() != 0:
                            FASTA_INCLUDED = FINAL_LOC[FINAL_LOC['PERCENTAGE'] >= Localization_identity_var.get()].reset_index(drop=True).drop('FASTA_UNIPROT', axis=1).drop_duplicates()

                            SUBCELLULAR_LOCALIZATION = pd.DataFrame({
                                "FASTA ID" : FASTA_INCLUDED["FASTA ID"],
                                "LOC" : FASTA_INCLUDED["LOC"],
                                "PERCENTAGE" : FASTA_INCLUDED["PERCENTAGE"],
                                "FASTA_SUBMITTED" : FASTA_INCLUDED["FASTA_SUBMITTED"]
                            }).drop_duplicates().reset_index(drop=True)

                            if len(FASTA_INCLUDED) == 0:
                                messagebox.showerror("Error", "No Relevant Data found")
                            else:
                                with open(f"{SAVE_FILE_Var.get()}/FIRST_LOCALIZATION.faa", "w") as faa_file:
                                    for FAA_PROTEIN, FAA_LOCATION, FAA_ID, FAA_PERC in zip(
                                            SUBCELLULAR_LOCALIZATION['FASTA_SUBMITTED'],
                                            SUBCELLULAR_LOCALIZATION['LOC'],
                                            SUBCELLULAR_LOCALIZATION['FASTA ID'],
                                            SUBCELLULAR_LOCALIZATION['PERCENTAGE']):
                                        faa_file.write(
                                            f">{FAA_ID} {FAA_LOCATION} {FAA_PERC}%" + "\n" + str(FAA_PROTEIN) + "\n")
                                    faa_file.close()
                                Localization_Proteins.config(
                                    text=f"{len(SUBCELLULAR_LOCALIZATION)}\nProtein(s)",
                                    fg="Red")
                                # messagebox.showinfo("Completed", "Localizaton Completed")
                                sleep(2)

                if Non_Homologs_var.get() == "True":

                    HUMAN_TITLE, HUMAN_FASTA, HUMAN_DES = read_submitted_fasta_file("./Data/Human_Reference_Genome.fasta")

                    FINAL_DATA = pd.DataFrame({

                    })

                    COUNTED_VAL_LIST = []
                    COUNTED_VAL_SEQ = []
                    COUNTED_VAL_ID = []

                    for single in data_list:
                        for single_fasta, single_title in zip(single["FASTA"], single["TITLE"]):
                            FINAL_DATA["HUMAN"] = HUMAN_FASTA
                            FINAL_DATA["SUBMITTED_SEQUENCE"] = single_fasta
                            FINAL_DATA["SUBMITTED ID"] = single_title
                            FINAL_DATA["LEN_HUMAN"] = FINAL_DATA["HUMAN"].apply(len)
                            FINAL_DATA["LEN_SUBMITTED"] = FINAL_DATA["SUBMITTED_SEQUENCE"].apply(len)
                            FILTERED_LENGTH = FINAL_DATA[FINAL_DATA["LEN_SUBMITTED"] <= FINAL_DATA["LEN_HUMAN"]].reset_index(drop=True)
                            FILTERED_LENGTH["RESULT"] = FILTERED_LENGTH.swifter.apply(lambda x: compute_matched([x["HUMAN"], x["SUBMITTED_SEQUENCE"]]), axis=1)
                            CLEANED = FILTERED_LENGTH[FILTERED_LENGTH["RESULT"] != "None"]
                            if Non_Homologs_identity_var.get() != 0:
                                IDENTITY_FILTERED_HOMOLOGS = CLEANED[CLEANED["RESULT"] <= Non_Homologs_identity_var.get()]
                                counts = IDENTITY_FILTERED_HOMOLOGS.value_counts(["SUBMITTED_SEQUENCE", "SUBMITTED ID"]).to_dict()
                                for COUNTED_VAL, COUNTED_NAMES in zip(counts.values(), counts.keys()):
                                    if Non_Homologs_perc_var.get() != 0:
                                        percentage = int(COUNTED_VAL / len(HUMAN_FASTA) * 100)
                                        if percentage >= Non_Homologs_perc_var.get():
                                            COUNTED_VAL_LIST.append(f"{percentage}%")
                                            COUNTED_VAL_SEQ.append(COUNTED_NAMES[0])
                                            COUNTED_VAL_ID.append(COUNTED_NAMES[1])

                    HOMOLOGS_DATA = pd.DataFrame({
                        "TITLE": COUNTED_VAL_ID,
                        "FASTA": COUNTED_VAL_SEQ,
                        "DES": COUNTED_VAL_LIST
                    }).drop_duplicates().reset_index(drop=True)

                    print(len(HOMOLOGS_DATA))

                    if len(HOMOLOGS_DATA) == 0:
                        messagebox.showerror("Error", "No Relevant Data found")
                    else:
                        with open(f"{SAVE_FILE_Var.get()}//NON_HOMOLOGS.faa", "w") as faa_file:
                            for FAA_ID,FAA_PROTEIN, FAA_PERC in zip(HOMOLOGS_DATA["TITLE"], HOMOLOGS_DATA["FASTA"],
                                                                     HOMOLOGS_DATA["DES"]):
                                faa_file.write(
                                    f">{FAA_ID} {FAA_PERC} NON HOMOLOGY" + "\n" + str(FAA_PROTEIN) + "\n")
                            faa_file.close()
                            Non_Homologs_Proteins.config(text=f"{len(HOMOLOGS_DATA)}\nProtein(s)", fg="Red")
                            sleep(2)

                if VFDB_var.get() == "True":

                    VFDB = pd.read_csv("./Data/VFDB.csv")

                    DATA_FRAME = pd.DataFrame({})
                    LIST_LOC = []

                    for single in data_list:
                        for single_fasta, single_title in zip(single["FASTA"], single["TITLE"]):
                            DATA_FRAME["VFDB_SEQUENCE"] = VFDB["FASTA"]
                            DATA_FRAME['SUBMITTED ID'] = single_title
                            DATA_FRAME["SUBMITTED_SEQUENCE"] = single_fasta
                            DATA_FRAME["LEN_VFDB"] = DATA_FRAME["VFDB_SEQUENCE"].apply(len)
                            DATA_FRAME["LEN_SUBMITTED"] = DATA_FRAME["SUBMITTED_SEQUENCE"].apply(len)
                            FILTERED_LENGTH = DATA_FRAME[
                                DATA_FRAME["LEN_SUBMITTED"] <= DATA_FRAME["LEN_VFDB"]].reset_index(drop=True)
                            try:
                                FILTERED_LENGTH["PERCENTAGE"] = FILTERED_LENGTH.apply(
                                    lambda x: compute_matched([x["VFDB_SEQUENCE"], x["SUBMITTED_SEQUENCE"]]),
                                    axis=1)
                                FINAL = FILTERED_LENGTH[
                                    FILTERED_LENGTH["PERCENTAGE"] == FILTERED_LENGTH['PERCENTAGE'].max()]
                                FINAL_SEP = FINAL[FINAL["LEN_VFDB"] == FINAL[["LEN_VFDB", "LEN_SUBMITTED"]][
                                    "LEN_VFDB"].min()].reset_index(drop=True)
                                LIST_LOC.append(FINAL_SEP)
                            except ValueError:
                                pass

                    FINAL_LOC = pd.concat(LIST_LOC)

                    FINAL_VFDB = FINAL_LOC[FINAL_LOC["PERCENTAGE"] >= VFDB_identity_var.get()].reset_index(drop=True)

                    VFDB_DATAFRAME = pd.DataFrame({
                        "SUBMITTED_SEQUENCE" : FINAL_VFDB["SUBMITTED_SEQUENCE"],
                        "SUBMITTED ID" : FINAL_VFDB["SUBMITTED ID"],
                        "PERCENTAGE" : FINAL_VFDB["PERCENTAGE"]
                    }).drop_duplicates().reset_index(drop=True)

                    if len(DATA_FRAME) == 0:
                        messagebox.showerror("Error", "No Relevant Data found")
                    else:
                        with open(f"{SAVE_FILE_Var.get()}//VFDB.faa", "w") as faa_file:
                            for FAA_PROTEIN, FAA_ID, FAA_PERC in zip(
                                    VFDB_DATAFRAME['SUBMITTED_SEQUENCE'],
                                    VFDB_DATAFRAME['SUBMITTED ID'],
                                    VFDB_DATAFRAME['PERCENTAGE']):
                                faa_file.write(f">{FAA_ID} {FAA_PERC}% Virulence" + "\n" + str(FAA_PROTEIN) + "\n")
                            faa_file.close()
                        VFDB_Proteins.config(
                            text=f"{len(VFDB_DATAFRAME.drop_duplicates())}\nProtein(s)",
                            fg="Red")
                        sleep(2)

                if Epitope_var.get() == "True":

                    EPI, EPI_FASTAS, EPI_ORGANISM, = read_submitted_fasta_file("./Data/Linear_Epitopes.faa")
                    CD8 = pd.read_csv("./Data/CD8_EPITOPES.csv").drop("Unnamed: 0", axis=1)
                    CD4 = pd.read_csv("./Data/CD4_EPITOPES.csv").drop("Unnamed: 0", axis=1)

                    EPI_DATA_FRAME = pd.DataFrame({"EPI": EPI_FASTAS})
                    ALL_EPITOPES_B_CELL = []

                    CD8_DATA_FRAME = pd.DataFrame({})
                    CD4_DATA_FRAME = pd.DataFrame({})

                    CD8_EPITOPES = []
                    CD4_EPITOPES = []

                    if Epitope_window_B_cell_var.get() == "ALL":
                        WINDOW_LENGTH = [12, 13, 14, 15, 16, 17, 18, 19, 20]
                    else:
                        WINDOW_LENGTH = [int(Epitope_window_B_cell_var.get())]

                    if Epitope_window_CD4_var.get() == "ALL":
                        CD4_WINDOW_LENGTH = [13, 14, 15, 16, 17, 18, 19, 20]
                    else:
                        CD4_WINDOW_LENGTH = [int(Epitope_window_CD4_var.get())]

                    if Epitope_window_CD8_var.get() == "ALL":
                        CD8_WINDOW_LENGTH = [8, 9, 10, 11, 12, 13]
                    else:
                        CD8_WINDOW_LENGTH = [int(Epitope_window_CD8_var.get())]

                    for single in data_list:
                        VIRULENCE_DATA_FRAME = pd.DataFrame({
                            "ID": single["TITLE"],
                            "FASTA": single["FASTA"],
                        }).drop_duplicates()

                        for VFDB_FASTA, IDs_B_CELL in zip(VIRULENCE_DATA_FRAME["FASTA"], VIRULENCE_DATA_FRAME["ID"]):
                            for WINDOW_LENGTHs in WINDOW_LENGTH:
                                FASTA_CHUNKS = [VFDB_FASTA[i: j] for i in range(len(VFDB_FASTA)) for j in
                                                range(i + 1, len(VFDB_FASTA) + 1) if
                                                len(VFDB_FASTA[i:j]) == WINDOW_LENGTHs]
                                for SINGLE_CHUNK in FASTA_CHUNKS:
                                    FILTERED_LEN = EPI_DATA_FRAME[
                                        EPI_DATA_FRAME["EPI"].apply(len) == WINDOW_LENGTHs].reset_index(drop=True)
                                    FILTERED_LEN["VIRULENCE CHUNK"] = SINGLE_CHUNK
                                    FILTERED_LEN["ID"] = IDs_B_CELL
                                    FILTERED_LEN["TYPE"] = "B CELL"
                                    FILTERED_LEN["RESULT"] = FILTERED_LEN.apply(
                                        lambda x: compute_matched([x["EPI"], x["VIRULENCE CHUNK"]]), axis=1)
                                    B_CELL_RESULT = FILTERED_LEN[FILTERED_LEN["RESULT"] >= Epitope_identity_B_cell_var.get()]
                                    ALL_EPITOPES_B_CELL.append(B_CELL_RESULT)

                    CONCAT_B_CELL = pd.concat(ALL_EPITOPES_B_CELL).reset_index(drop=True)

                    for B_CELL_EPITOPE, IDs_T_CELL in zip(VIRULENCE_DATA_FRAME["FASTA"],
                                                          VIRULENCE_DATA_FRAME["ID"]):

                        for CD8_WINDOW_LENGTHs in CD8_WINDOW_LENGTH:
                            FASTA_CHUNKS_CD8 = [B_CELL_EPITOPE[i: j] for i in range(len(B_CELL_EPITOPE)) for j in
                                                range(i + 1, len(B_CELL_EPITOPE) + 1) if
                                                len(B_CELL_EPITOPE[i:j]) == CD8_WINDOW_LENGTHs]

                        for CD4_WINDOW_LENGTHs in CD4_WINDOW_LENGTH:
                            FASTA_CHUNKS_CD4 = [B_CELL_EPITOPE[i: j] for i in range(len(B_CELL_EPITOPE)) for j in
                                                range(i + 1, len(B_CELL_EPITOPE) + 1) if
                                                len(B_CELL_EPITOPE[i:j]) == CD4_WINDOW_LENGTHs]

                            for SINGLE_CHUNK_CD8 in FASTA_CHUNKS_CD8:
                                CD8_DATA_FRAME["ID"] = IDs_T_CELL
                                CD8_DATA_FRAME["TYPE"] = "CD8"
                                CD8_DATA_FRAME["ALLEL"] = CD8["ALLEL"]
                                CD8_DATA_FRAME["EPI"] = CD8["FASTA"]
                                FILTERED_LEN_CD8 = CD8_DATA_FRAME[
                                    CD8_DATA_FRAME["EPI"].apply(len) == CD8_WINDOW_LENGTHs].reset_index(drop=True)
                                CD8_DATA_FRAME["VIRULENCE CHUNK"] = SINGLE_CHUNK_CD8
                                CD8_DATA_FRAME["RESULT"] = CD8_DATA_FRAME.apply(
                                    lambda x: compute_matched([x["EPI"], x["VIRULENCE CHUNK"]]), axis=1)
                                CD8_RESULT = CD8_DATA_FRAME[CD8_DATA_FRAME["RESULT"] >= Epitope_identity_CD8_var.get()]
                                CD8_EPITOPES.append(CD8_RESULT)

                            for SINGLE_CHUNK_CD4 in FASTA_CHUNKS_CD4:
                                CD4_DATA_FRAME["ID"] = IDs_T_CELL
                                CD4_DATA_FRAME["TYPE"] = "CD4"
                                CD4_DATA_FRAME["ALLEL"] = CD4["ALLEL"]
                                CD4_DATA_FRAME["EPI"] = CD4["FASTA"]
                                FILTERED_LEN_CD4 = CD4_DATA_FRAME[
                                    CD4_DATA_FRAME["EPI"].apply(len) == CD4_WINDOW_LENGTHs].reset_index(drop=True)
                                CD4_DATA_FRAME["VIRULENCE CHUNK"] = SINGLE_CHUNK_CD4
                                CD4_DATA_FRAME["RESULT"] = CD4_DATA_FRAME.apply(
                                    lambda x: compute_matched([x["EPI"], x["VIRULENCE CHUNK"]]), axis=1)
                                CD4_RESULT = CD4_DATA_FRAME[CD4_DATA_FRAME["RESULT"] >= Epitope_identity_CD4_var.get()]
                                CD4_EPITOPES.append(CD4_RESULT)

                    CD8_EPITOPES_DATAFRAME = pd.concat(CD8_EPITOPES).reset_index(drop=True).dropna()
                    CD4_EPITOPES_DATAFRAME = pd.concat(CD4_EPITOPES).reset_index(drop=True).dropna()

                    T_CELL_FILTERED = pd.concat([CD8_EPITOPES_DATAFRAME, CD4_EPITOPES_DATAFRAME])

                    T_CELL_DATAFRAME = pd.DataFrame({
                        "EPI" : T_CELL_FILTERED["EPI"],
                        "ID" : T_CELL_FILTERED["ID"],
                        "RESULT" : T_CELL_FILTERED["RESULT"],
                        "TYPE" : T_CELL_FILTERED["TYPE"],
                        "ALLEL" : T_CELL_FILTERED["ALLEL"]
                    }).drop_duplicates().reset_index(drop=True)

                    B_CELL_DATAFRAME = pd.DataFrame({
                        "EPI" : CONCAT_B_CELL["EPI"],
                        "ID" : CONCAT_B_CELL["ID"],
                        "RESULT" : CONCAT_B_CELL["RESULT"],
                        "TYPE" : CONCAT_B_CELL["TYPE"]
                    }).drop_duplicates().reset_index(drop=True)

                    with open(f"{SAVE_FILE_Var.get()}//T_EPITOPES.faa", "w") as faa_file_first:
                        for FAA_PROTEIN, FAA_ID, FAA_PERC, FAA_NAME, FAA_TITLE in zip(
                                T_CELL_DATAFRAME['EPI'],
                                T_CELL_DATAFRAME['ID'],
                                T_CELL_DATAFRAME['RESULT'],
                                T_CELL_DATAFRAME['TYPE'],
                                T_CELL_DATAFRAME["ALLEL"]):
                            faa_file_first.write(
                                f">{FAA_ID} {int(FAA_PERC)}% {FAA_NAME} {FAA_TITLE}" + "\n" + str(FAA_PROTEIN) + "\n")
                        faa_file_first.close()

                    with open(f"{SAVE_FILE_Var.get()}//B_EPITOPES.faa", "w") as faa_file_second:
                        for FAA_PROTEIN, FAA_ID, FAA_PERC, FAA_NAME in zip(
                                B_CELL_DATAFRAME['EPI'],
                                B_CELL_DATAFRAME['ID'],
                                B_CELL_DATAFRAME['RESULT'],
                                B_CELL_DATAFRAME['TYPE']):
                            faa_file_second.write(
                                f">{FAA_ID} {int(FAA_PERC)}% {FAA_NAME} {FAA_TITLE}" + "\n" + str(FAA_PROTEIN) + "\n")
                        faa_file_second.close()

                    TOTAL_EPI = len(T_CELL_DATAFRAME) + len(B_CELL_DATAFRAME)

                    Epitope_peptides.config(text=f"{TOTAL_EPI}\nEpitope(s)", fg="Red")
                    sleep(2)

                    Localization_Proteins.config(text=f"0\nProtein(s)", fg=Label_Color)
                    Non_Homologs_Proteins.config(text=f"0\nProtein(s)", fg=Label_Color)
                    VFDB_Proteins.config(text=f"0\nProtein(s)", fg=Label_Color)
                    Epitope_peptides.config(text=f"0\nEpitope(s)", fg=Label_Color)

                Update_Progress(Bar_Value[0] * len(Progress_bar_list))

                sleep(2)
                Localization_Proteins.config(text=f"0\nProtein(s)",fg=Label_Color)
                Non_Homologs_Proteins.config(text=f"0\nProtein(s)",fg=Label_Color)
                VFDB_Proteins.config(text=f"0\nProtein(s)",fg=Label_Color)
                Epitope_peptides.config(text=f"0\nEpitope(s)",fg=Label_Color)

                Update_Progress(0)

            else:
                messagebox.showwarning("Warning", "No Folder Selected To Place Data")

        else:

            Job_Status.config(text="Job Status: Running", fg="Red")

            if SAVE_FILE_Var.get() != "":

                try:

                    Update_Progress(5)

                    LOCATION_LIST = []
                    BIT_SCORE_LIST = []
                    PROTEIN_ID_LIST = []
                    FASTA_LIST = []

                    TITLE, FASTA, DESCRIPTION = read_submitted_fasta_file(f"./Data/FINAL_FILES/{FILE_NAME_VAR.get()}")

                    for TITLE_READED, FASTA_READED, DESCRIPTION_READED in zip(TITLE, FASTA, DESCRIPTION):
                        info = " ".join(DESCRIPTION_READED.split(" ")[1:-1])
                        score = DESCRIPTION_READED.split(" ")[-1]

                        LOCATION_LIST.append(info)
                        BIT_SCORE_LIST.append(int(score))
                        PROTEIN_ID_LIST.append(TITLE_READED)
                        FASTA_LIST.append(FASTA_READED)

                    LOCATION_DATA = pd.DataFrame({
                        "PROTEIN ID": PROTEIN_ID_LIST,
                        "BIT SCORE": BIT_SCORE_LIST,
                        "LOCATION": LOCATION_LIST,
                        "FASTA": FASTA_LIST
                    }).reset_index(drop=True)

                    if Reliability_var.get() != 0:
                        mean_filtered = LOCATION_DATA[LOCATION_DATA["BIT SCORE"] > Reliability_var.get()]
                        remove_list = ["secreted", "outer membrane", "fimbrium"]
                        NEW_DATA_FRAME = mean_filtered[mean_filtered["LOCATION"].isin(remove_list)].reset_index(
                            drop=True)
                        REALIABILITY_VAR.set(Reliability_var.get())

                        DATA_FINAL_LOC = pd.DataFrame({})
                        LIST_LOC = []

                        for single_Localization in data_list:
                            for single_Localization_fasta, single_Localization_id in zip(single_Localization["FASTA"],
                                                                                         single_Localization["TITLE"]):
                                DATA_FINAL_LOC["LOC"] = NEW_DATA_FRAME["LOCATION"]
                                DATA_FINAL_LOC["FASTA_UNIPROT"] = NEW_DATA_FRAME["FASTA"]
                                DATA_FINAL_LOC['FASTA_SUBMITTED'] = single_Localization_fasta
                                DATA_FINAL_LOC['FASTA ID'] = single_Localization_id
                                DATA_FINAL_LOC["LEN_UNIPROT"] = DATA_FINAL_LOC["FASTA_UNIPROT"].apply(len)
                                DATA_FINAL_LOC["LEN_SUBMITTED"] = DATA_FINAL_LOC['FASTA_SUBMITTED'].apply(len)
                                FILTERED_LENGTH = DATA_FINAL_LOC[
                                    DATA_FINAL_LOC["LEN_SUBMITTED"] <= DATA_FINAL_LOC["LEN_UNIPROT"]].reset_index(drop=True)
                                try:
                                    FILTERED_LENGTH['PERCENTAGE'] = FILTERED_LENGTH.swifter.apply(
                                        lambda x: compute_matched([x["FASTA_UNIPROT"], x['FASTA_SUBMITTED']]), axis=1)
                                    FINAL = FILTERED_LENGTH[
                                        FILTERED_LENGTH['PERCENTAGE'] == FILTERED_LENGTH['PERCENTAGE'].max()]
                                    FINAL_SEP = FINAL[FINAL["LEN_UNIPROT"] == FINAL[["LEN_UNIPROT", "LEN_SUBMITTED"]][
                                        "LEN_UNIPROT"].min()].reset_index(drop=True)
                                    LIST_LOC.append(FINAL_SEP)
                                except ValueError:
                                    pass
                        FINAL_LOC = pd.concat(LIST_LOC)

                        if Localization_identity_var.get() != 0:
                            FASTA_INCLUDED = FINAL_LOC[
                                FINAL_LOC['PERCENTAGE'] >= Localization_identity_var.get()].reset_index(drop=True)

                            SUBCELLULAR_LOCALIZATION = pd.DataFrame({
                                "FASTA ID" : FASTA_INCLUDED["FASTA ID"],
                                "LOC" : FASTA_INCLUDED["LOC"],
                                "PERCENTAGE" : FASTA_INCLUDED["PERCENTAGE"],
                                "FASTA_SUBMITTED" : FASTA_INCLUDED["FASTA_SUBMITTED"]
                            }).drop_duplicates().reset_index(drop=True)

                            if len(FASTA_INCLUDED) == 0:
                                messagebox.showerror("Error", "No Relevant Data found")
                            else:
                                with open(f"{SAVE_FILE_Var.get()}/FIRST_LOCALIZATION_{len(SUBCELLULAR_LOCALIZATION)}.faa", "w") as faa_file:
                                    for FAA_PROTEIN, FAA_LOCATION, FAA_ID, FAA_PERC in zip(
                                            SUBCELLULAR_LOCALIZATION['FASTA_SUBMITTED'],
                                            SUBCELLULAR_LOCALIZATION['LOC'],
                                            SUBCELLULAR_LOCALIZATION['FASTA ID'],
                                            SUBCELLULAR_LOCALIZATION['PERCENTAGE']):
                                        faa_file.write(
                                            f">{FAA_ID} {FAA_LOCATION} {FAA_PERC}%" + "\n" + str(FAA_PROTEIN) + "\n")
                                    faa_file.close()
                                Localization_Proteins.config(
                                    text=f"{len(SUBCELLULAR_LOCALIZATION)}\nProtein(s)",
                                    fg="Red")
                                # messagebox.showinfo("Completed", "Localizaton Completed")
                                sleep(2)

                    Update_Progress(25)

                    HUMAN_TITLE, HUMAN_FASTA, HUMAN_DES = read_submitted_fasta_file(
                        "./Data/Human_Reference_Genome.fasta")

                    FINAL_DATA = pd.DataFrame({

                    })

                    COUNTED_VAL_LIST = []
                    COUNTED_VAL_SEQ = []
                    COUNTED_VAL_ID = []

                    for single_fasta, single_title in zip(FASTA_INCLUDED['FASTA_SUBMITTED'],FASTA_INCLUDED['FASTA ID']):
                        FINAL_DATA["HUMAN"] = HUMAN_FASTA
                        FINAL_DATA["SUBMITTED_SEQUENCE"] = single_fasta
                        FINAL_DATA["SUBMITTED ID"] = single_title
                        FINAL_DATA["LEN_HUMAN"] = FINAL_DATA["HUMAN"].apply(len)
                        FINAL_DATA["LEN_SUBMITTED"] = FINAL_DATA["SUBMITTED_SEQUENCE"].apply(len)
                        FILTERED_LENGTH = FINAL_DATA[
                            FINAL_DATA["LEN_SUBMITTED"] <= FINAL_DATA["LEN_HUMAN"]].reset_index(drop=True)
                        FILTERED_LENGTH["RESULT"] = FILTERED_LENGTH.swifter.apply(
                            lambda x: compute_matched([x["HUMAN"], x["SUBMITTED_SEQUENCE"]]), axis=1)
                        CLEANED = FILTERED_LENGTH[FILTERED_LENGTH["RESULT"] != "None"]
                        if Non_Homologs_identity_var.get() != 0:
                            IDENTITY_FILTERED_HOMOLOGS = CLEANED[
                                CLEANED["RESULT"] <= Non_Homologs_identity_var.get()]
                            counts = IDENTITY_FILTERED_HOMOLOGS.value_counts(
                                ["SUBMITTED_SEQUENCE", "SUBMITTED ID"]).to_dict()
                            for COUNTED_VAL, COUNTED_NAMES in zip(counts.values(), counts.keys()):
                                if Non_Homologs_perc_var.get() != 0:
                                    percentage = int(COUNTED_VAL / len(HUMAN_FASTA) * 100)
                                    if percentage >= Non_Homologs_perc_var.get():
                                        COUNTED_VAL_LIST.append(f"{percentage}%")
                                        COUNTED_VAL_SEQ.append(COUNTED_NAMES[0])
                                        COUNTED_VAL_ID.append(COUNTED_NAMES[1])

                    HOMOLOGS_DATA = pd.DataFrame({
                        "TITLE": COUNTED_VAL_ID,
                        "FASTA": COUNTED_VAL_SEQ,
                        "DES": COUNTED_VAL_LIST
                    }).drop_duplicates().reset_index(drop=True)

                    if len(FINAL_DATA) == 0:
                        messagebox.showerror("Error", "No Relevant Data found")
                    else:
                        with open(f"{SAVE_FILE_Var.get()}//NON_HOMOLOGS_{len(HOMOLOGS_DATA)}.faa", "w") as faa_file:
                            for FAA_ID,FAA_PROTEIN, FAA_PERC in zip(HOMOLOGS_DATA["TITLE"], HOMOLOGS_DATA["FASTA"],
                                                                     HOMOLOGS_DATA["DES"]):
                                faa_file.write(
                                    f">{FAA_ID} {FAA_PERC} NON HOMOLOGY" + "\n" + str(FAA_PROTEIN) + "\n")
                            faa_file.close()

                    Non_Homologs_Proteins.config(text=f"{len(HOMOLOGS_DATA)}\nProtein(s)", fg="Red")
                    sleep(2)

                    Update_Progress(50)

                    VFDB = pd.read_csv("./Data/VFDB.csv")

                    DATA_FRAME = pd.DataFrame({})
                    LIST_VFDB = []

                    for single_fasta, single_title in zip(COUNTED_VAL_SEQ, COUNTED_VAL_ID):
                        DATA_FRAME["VFDB_SEQUENCE"] = VFDB["FASTA"]
                        DATA_FRAME['SUBMITTED ID'] = single_title
                        DATA_FRAME["SUBMITTED_SEQUENCE"] = single_fasta
                        DATA_FRAME["LEN_VFDB"] = DATA_FRAME["VFDB_SEQUENCE"].apply(len)
                        DATA_FRAME["LEN_SUBMITTED"] = DATA_FRAME["SUBMITTED_SEQUENCE"].apply(len)
                        FILTERED_LENGTH = DATA_FRAME[
                            DATA_FRAME["LEN_SUBMITTED"] <= DATA_FRAME["LEN_VFDB"]].reset_index(drop=True)
                        try:
                            FILTERED_LENGTH["PERCENTAGE"] = FILTERED_LENGTH.apply(
                                lambda x: compute_matched([x["VFDB_SEQUENCE"], x["SUBMITTED_SEQUENCE"]]),
                                axis=1)
                            FINAL = FILTERED_LENGTH[
                                FILTERED_LENGTH["PERCENTAGE"] == FILTERED_LENGTH['PERCENTAGE'].max()]
                            FINAL_SEP = FINAL[FINAL["LEN_VFDB"] == FINAL[["LEN_VFDB", "LEN_SUBMITTED"]][
                                "LEN_VFDB"].min()].reset_index(drop=True)
                            LIST_VFDB.append(FINAL_SEP)
                        except ValueError:
                            pass

                    FINAL_VIRULENCE = pd.concat(LIST_VFDB)

                    FINAL_VFDB = FINAL_VIRULENCE[FINAL_VIRULENCE["PERCENTAGE"] >= VFDB_identity_var.get()].reset_index(drop=True)

                    VFDB_DATAFRAME = pd.DataFrame({
                        "SUBMITTED_SEQUENCE": FINAL_VFDB["SUBMITTED_SEQUENCE"],
                        "SUBMITTED ID": FINAL_VFDB["SUBMITTED ID"],
                        "PERCENTAGE": FINAL_VFDB["PERCENTAGE"]
                    }).drop_duplicates().reset_index(drop=True)

                    if len(DATA_FRAME) == 0:
                        messagebox.showerror("Error", "No Relevant Data found")
                    else:
                        with open(f"{SAVE_FILE_Var.get()}//VFDB_{len(VFDB_DATAFRAME)}.faa", "w") as faa_file:
                            for FAA_PROTEIN, FAA_ID, FAA_PERC in zip(
                                    VFDB_DATAFRAME['SUBMITTED_SEQUENCE'],
                                    VFDB_DATAFRAME['SUBMITTED ID'],
                                    VFDB_DATAFRAME['PERCENTAGE']):
                                faa_file.write(f">{FAA_ID} {FAA_PERC}% Virulence" + "\n" + str(FAA_PROTEIN) + "\n")
                            faa_file.close()

                        with open(f"{SAVE_FILE_Var.get()}//PVCs_{len(VFDB_DATAFRAME)}.faa", "w") as faa_file:
                            for FAA_PROTEIN, FAA_ID, FAA_PERC in zip(
                                    VFDB_DATAFRAME['SUBMITTED_SEQUENCE'],
                                    VFDB_DATAFRAME['SUBMITTED ID'],
                                    VFDB_DATAFRAME['PERCENTAGE']):
                                faa_file.write(f">{FAA_ID} {FAA_PERC}% Virulence" + "\n" + str(FAA_PROTEIN) + "\n")
                            faa_file.close()

                        VFDB_Proteins.config(
                            text=f"{len(VFDB_DATAFRAME.drop_duplicates())}\nProtein(s)",
                            fg="Red")

                        PVCs.config(
                            text=f"{len(VFDB_DATAFRAME.drop_duplicates())}\nPVC(s)",
                            fg="Green")

                        sleep(2)

                    Update_Progress(75)

                    EPI, EPI_FASTAS, EPI_ORGANISM, = read_submitted_fasta_file("./Data/Linear_Epitopes.faa")
                    CD8 = pd.read_csv("./Data/CD8_EPITOPES.csv").drop("Unnamed: 0", axis=1)
                    CD4 = pd.read_csv("./Data/CD4_EPITOPES.csv").drop("Unnamed: 0", axis=1)

                    EPI_DATA_FRAME = pd.DataFrame({"EPI": EPI_FASTAS})
                    ALL_EPITOPES_B_CELL = []

                    CD8_DATA_FRAME = pd.DataFrame({})
                    CD4_DATA_FRAME = pd.DataFrame({})

                    CD8_EPITOPES = []
                    CD4_EPITOPES = []

                    VIRULENCE_DATA_FRAME = pd.DataFrame({
                        "ID": FINAL_VFDB['SUBMITTED ID'],
                        "FASTA": FINAL_VFDB['SUBMITTED_SEQUENCE'],
                    }).drop_duplicates()

                    if Epitope_window_B_cell_var.get() == "ALL":
                        WINDOW_LENGTH = [12, 13, 14, 15, 16, 17, 18, 19, 20]
                    else:
                        WINDOW_LENGTH = [int(Epitope_window_B_cell_var.get())]

                    if Epitope_window_CD4_var.get() == "ALL":
                        CD4_WINDOW_LENGTH = [13, 14, 15, 16, 17, 18, 19, 20]
                    else:
                        CD4_WINDOW_LENGTH = [int(Epitope_window_CD4_var.get())]

                    if Epitope_window_CD8_var.get() == "ALL":
                        CD8_WINDOW_LENGTH = [8, 9, 10, 11, 12, 13]
                    else:
                        CD8_WINDOW_LENGTH = [int(Epitope_window_CD8_var.get())]

                    for WINDOW_LENGTHs in WINDOW_LENGTH:
                        for VFDB_FASTA, IDs_B_CELL in zip(VIRULENCE_DATA_FRAME["FASTA"], VIRULENCE_DATA_FRAME["ID"]):
                            FASTA_CHUNKS = [VFDB_FASTA[i: j] for i in range(len(VFDB_FASTA)) for j in
                                            range(i + 1, len(VFDB_FASTA) + 1) if
                                            len(VFDB_FASTA[i:j]) == WINDOW_LENGTHs]
                        for SINGLE_CHUNK in FASTA_CHUNKS:
                            FILTERED_LEN = EPI_DATA_FRAME[
                                EPI_DATA_FRAME["EPI"].apply(len) == WINDOW_LENGTHs].reset_index(drop=True)
                            FILTERED_LEN["VIRULENCE CHUNK"] = SINGLE_CHUNK
                            FILTERED_LEN["ID"] = IDs_B_CELL
                            FILTERED_LEN["TYPE"] = "B CELL"
                            FILTERED_LEN["RESULT"] = FILTERED_LEN.apply(
                                lambda x: compute_matched([x["EPI"], x["VIRULENCE CHUNK"]]), axis=1)
                            B_CELL_RESULT = FILTERED_LEN[FILTERED_LEN["RESULT"] >= Epitope_identity_B_cell_var.get()]
                            ALL_EPITOPES_B_CELL.append(B_CELL_RESULT)

                        CONCAT_B_CELL = pd.concat(ALL_EPITOPES_B_CELL).reset_index(drop=True)

                    for B_CELL_EPITOPE, IDs_T_CELL in zip(VIRULENCE_DATA_FRAME["FASTA"],
                                                          VIRULENCE_DATA_FRAME["ID"]):

                        for CD8_WINDOW_LENGTHs in CD8_WINDOW_LENGTH:
                            FASTA_CHUNKS_CD8 = [B_CELL_EPITOPE[i: j] for i in range(len(B_CELL_EPITOPE)) for j in
                                                range(i + 1, len(B_CELL_EPITOPE) + 1) if
                                                len(B_CELL_EPITOPE[i:j]) == CD8_WINDOW_LENGTHs]

                        for CD4_WINDOW_LENGTHs in CD4_WINDOW_LENGTH:
                            FASTA_CHUNKS_CD4 = [B_CELL_EPITOPE[i: j] for i in range(len(B_CELL_EPITOPE)) for j in
                                                range(i + 1, len(B_CELL_EPITOPE) + 1) if
                                                len(B_CELL_EPITOPE[i:j]) == CD4_WINDOW_LENGTHs]

                            for SINGLE_CHUNK_CD8 in FASTA_CHUNKS_CD8:
                                CD8_DATA_FRAME["ID"] = IDs_T_CELL
                                CD8_DATA_FRAME["TYPE"] = "CD8"
                                CD8_DATA_FRAME["ALLEL"] = CD8["ALLEL"]
                                CD8_DATA_FRAME["EPI"] = CD8["FASTA"]
                                FILTERED_LEN_CD8 = CD8_DATA_FRAME[
                                    CD8_DATA_FRAME["EPI"].apply(len) == CD8_WINDOW_LENGTHs].reset_index(drop=True)
                                CD8_DATA_FRAME["VIRULENCE CHUNK"] = SINGLE_CHUNK_CD8
                                CD8_DATA_FRAME["RESULT"] = CD8_DATA_FRAME.apply(
                                    lambda x: compute_matched([x["EPI"], x["VIRULENCE CHUNK"]]), axis=1)
                                CD8_RESULT = CD8_DATA_FRAME[CD8_DATA_FRAME["RESULT"] >= Epitope_identity_CD8_var.get()]
                                CD8_EPITOPES.append(CD8_RESULT)

                            for SINGLE_CHUNK_CD4 in FASTA_CHUNKS_CD4:
                                CD4_DATA_FRAME["ID"] = IDs_T_CELL
                                CD4_DATA_FRAME["TYPE"] = "CD4"
                                CD4_DATA_FRAME["ALLEL"] = CD4["ALLEL"]
                                CD4_DATA_FRAME["EPI"] = CD4["FASTA"]
                                FILTERED_LEN_CD4 = CD4_DATA_FRAME[
                                    CD4_DATA_FRAME["EPI"].apply(len) == CD4_WINDOW_LENGTHs].reset_index(drop=True)
                                CD4_DATA_FRAME["VIRULENCE CHUNK"] = SINGLE_CHUNK_CD4
                                CD4_DATA_FRAME["RESULT"] = CD4_DATA_FRAME.apply(
                                    lambda x: compute_matched([x["EPI"], x["VIRULENCE CHUNK"]]), axis=1)
                                CD4_RESULT = CD4_DATA_FRAME[CD4_DATA_FRAME["RESULT"] >= Epitope_identity_CD4_var.get()]
                                CD4_EPITOPES.append(CD4_RESULT)

                    CD8_EPITOPES_DATAFRAME = pd.concat(CD8_EPITOPES).reset_index(drop=True).dropna()
                    CD4_EPITOPES_DATAFRAME = pd.concat(CD4_EPITOPES).reset_index(drop=True).dropna()

                    T_CELL_FILTERED = pd.concat([CD8_EPITOPES_DATAFRAME, CD4_EPITOPES_DATAFRAME])

                    T_CELL_DATAFRAME = pd.DataFrame({
                        "EPI": T_CELL_FILTERED["EPI"],
                        "ID": T_CELL_FILTERED["ID"],
                        "RESULT": T_CELL_FILTERED["RESULT"],
                        "TYPE": T_CELL_FILTERED["TYPE"],
                        "ALLEL": T_CELL_FILTERED["ALLEL"]
                    }).drop_duplicates().reset_index(drop=True)

                    B_CELL_DATAFRAME = pd.DataFrame({
                        "EPI": CONCAT_B_CELL["EPI"],
                        "ID": CONCAT_B_CELL["ID"],
                        "RESULT": CONCAT_B_CELL["RESULT"],
                        "TYPE": CONCAT_B_CELL["TYPE"]
                    }).drop_duplicates().reset_index(drop=True)

                    with open(f"{SAVE_FILE_Var.get()}//T_EPITOPES_{len(T_CELL_DATAFRAME)}.faa", "w") as faa_file_first:
                        for FAA_PROTEIN, FAA_ID, FAA_PERC, FAA_NAME, FAA_TITLE in zip(
                                T_CELL_DATAFRAME['EPI'],
                                T_CELL_DATAFRAME['ID'],
                                T_CELL_DATAFRAME['RESULT'],
                                T_CELL_DATAFRAME['TYPE'],
                                T_CELL_DATAFRAME["ALLEL"]):
                            faa_file_first.write(
                                f">{FAA_ID} {int(FAA_PERC)}% {FAA_NAME} {FAA_TITLE}" + "\n" + str(FAA_PROTEIN) + "\n")
                        faa_file_first.close()

                    with open(f"{SAVE_FILE_Var.get()}//B_EPITOPES_{len(B_CELL_DATAFRAME)}.faa", "w") as faa_file_second:
                        for FAA_PROTEIN, FAA_ID, FAA_PERC, FAA_NAME in zip(
                                B_CELL_DATAFRAME['EPI'],
                                B_CELL_DATAFRAME['ID'],
                                B_CELL_DATAFRAME['RESULT'],
                                B_CELL_DATAFRAME['TYPE']):
                            faa_file_second.write(
                                f">{FAA_ID} {int(FAA_PERC)}% {FAA_NAME} {FAA_TITLE}" + "\n" + str(FAA_PROTEIN) + "\n")
                        faa_file_second.close()

                    TOTAL_EPI = len(T_CELL_DATAFRAME) + len(B_CELL_DATAFRAME)

                    Epitope_peptides.config(text=f"{TOTAL_EPI}\nEpitope(s)", fg="Red")
                    sleep(2)

                    Update_Progress(100)

                    sleep(2)

                    Update_Progress(0)

                    messagebox.showinfo("Completed", "Analysis Completed")

                    Localization_Proteins.config(text=f"0\nProtein(s)", fg=Label_Color)
                    Non_Homologs_Proteins.config(text=f"0\nProtein(s)", fg=Label_Color)
                    VFDB_Proteins.config(text=f"0\nProtein(s)", fg=Label_Color)
                    PVCs.config(text=f"0\nPVC(s)", fg=Label_Color)
                    Epitope_peptides.config(text=f"0\nEpitope(s)", fg=Label_Color)

                    sleep(2)

                except:
                    messagebox.showerror("Error", "Make Sure You Choose The Parameters Correctly")
                    Update_Progress(0)
                    Localization_Proteins.config(text=f"0\nProtein(s)", fg=Label_Color)
                    Non_Homologs_Proteins.config(text=f"0\nProtein(s)", fg=Label_Color)
                    VFDB_Proteins.config(text=f"0\nProtein(s)", fg=Label_Color)
                    PVCs.config(text=f"0\nPVC(s)", fg=Label_Color)
                    Epitope_peptides.config(text=f"0\nEpitope(s)", fg=Label_Color)

            else:
                messagebox.showwarning("Warning", "No Folder Selected To Place Data")

    Job_Status.config(text="Job Status: Not Running", fg=Label_Color)


Save_file_btn = Button(Submission_method_frame, image=SAVE_BUTTON, command=SAVE_FILE, bd=0)
Save_file_btn.place(relx=0.79, rely=0.72)

Progress_bar_Frame = Frame(window, width=width_frames, height=30, bg=Second_Color, bd=0)
Progress_bar_Frame.grid(row=1, sticky="w")

my_progress = ttk.Progressbar(Progress_bar_Frame, orient=HORIZONTAL, length=900, mode="determinate")
my_progress.place(relx=0.08, rely=0.1)

Must_Evalutate_Frame = Frame(window, width=width_frames, height=30, bg=Second_Color, bd=0)
Must_Evalutate_Frame.grid(row=2, sticky="w")

Must_Evalutate_var = StringVar()
Must_Evalutate = Checkbutton(Must_Evalutate_Frame, variable=Must_Evalutate_var,text="Must Evaluate" ,onvalue="True", offvalue="False", bg=Second_Color, activebackground=Second_Color, command=MUST_EVALUATE, font=(font_family, 11))
Must_Evalutate.deselect()
Must_Evalutate.place(relx=0.015, rely=0.1)

Must_Evalutate_Details = Label(Must_Evalutate_Frame, text="If must evaluate is checked do not uncheck any method",bg=Second_Color, fg="Red", font=(font_family, 10))
Must_Evalutate_Details.place(relx=0.15, rely=0.14)

Example_Btn = Button(Submission_method_frame, text="Example Sequence",bg=Second_Color, fg=Buttons_Colors, font=(font_family, 10), bd=0,background=First_Color, command=FILL_EXAMPLE)
Example_Btn.place(relx=0.675, rely=0.3)

Localization_Frame = Frame(window, width=width_frames, height=60, bg=Second_Color, bd=0)
Localization_Frame.grid(row=3, sticky="w")

Localization_var = StringVar()
Localization = Checkbutton(Localization_Frame, variable=Localization_var,text="Localization", onvalue="True", offvalue="False", bg=Second_Color, activebackground=Second_Color, font=(font_family, 11))
Localization.deselect()
Localization.place(relx=0.015, rely=0.38)

Localization_Proteins = Label(Localization_Frame, text="0\nProtein(s)", bg=Second_Color, fg=Label_Color, font=(font_family, 12))
Localization_Proteins.place(relx=0.87, rely=0.15)

Localization_file_name = StringVar()

def get_single_strain(single):
    if single[0][0].isupper() == True:
        word = single[0]+"_"+single[1]
        files_list = os.listdir("./Data/FINAL_FILES")
        for i in files_list:
            name = "_".join(i.split(".")[0].split("_")[0:2])
            if name == word:
                return i
    else:
        word = single[1]+"_"+single[0]
        files_list = os.listdir("./Data/FINAL_FILES")
        for i in files_list:
            name = "_".join(i.split(".")[0].split("_")[0:2])
            if name == word:
                return i

# get_names() function will fetch all of the names of the bacteria present in the directory "Bacteria"
def get_names():
    refined_list = []
    files_list = os.listdir("./Data/FINAL_FILES")
    for i in files_list:
        refined_list.append({
            "First_Name": i.split("_")[0],
            "Second_Name": i.split("_")[1].split(".")[0]
        })
    return refined_list

genus = get_names()

#short variable represents total number of bacteria genus.
shorted_list = []
for i in genus:
    shorted_list.append(i["First_Name"])
short = list(set(shorted_list))


LOC_TREE_var = StringVar()
FILE_NAME_VAR = StringVar()
REALIABILITY_VAR = IntVar()

def Localization_select(event):

    final_list = []
    genus_selected = Localization_combo_box_1_var.get()
    second_name_filtered = []

    for x in genus:
        if genus_selected == x["First_Name"]:
            second_name_filtered.append(x["Second_Name"])
            final_list.append(x["First_Name"])

    Localization_combo_box_2.config(values=second_name_filtered)

    for y in second_name_filtered:
        if y == Localization_combo_box_2.get():
            final_list.append(y)

    single = list(set(final_list))
    if len(single) == 2:
        file_name = get_single_strain(single)
        FILE_NAME_VAR.set(file_name)


        LOCATION_LIST = []
        BIT_SCORE_LIST = []
        PROTEIN_ID_LIST = []
        FASTA_LIST = []

        TITLE, FASTA, DESCRIPTION = read_submitted_fasta_file(f"./Data/FINAL_FILES/{file_name}")

        for TITLE_READED, FASTA_READED, DESCRIPTION_READED in zip(TITLE, FASTA, DESCRIPTION):
            info = " ".join(DESCRIPTION_READED.split(" ")[1:-1])
            score = DESCRIPTION_READED.split(" ")[-1]

            LOCATION_LIST.append(info)
            BIT_SCORE_LIST.append(int(score))
            PROTEIN_ID_LIST.append(TITLE_READED)
            FASTA_LIST.append(FASTA_READED)

        LOCATION_DATA = pd.DataFrame({
            "PROTEIN ID": PROTEIN_ID_LIST,
            "BIT SCORE": BIT_SCORE_LIST,
            "LOCATION": LOCATION_LIST,
            "FASTA": FASTA_LIST
        }).reset_index(drop=True)

        if Reliability_var.get() != 0:

            mean_filtered = LOCATION_DATA[LOCATION_DATA["BIT SCORE"] > Reliability_var.get()]
            remove_list = ["secreted", "outer membrane", "fimbrium"]
            NEW_DATA_FRAME = mean_filtered[mean_filtered["LOCATION"].isin(remove_list)].reset_index(drop=True)
            REALIABILITY_VAR.set(Reliability_var.get())

            if len(NEW_DATA_FRAME) == 0:
                messagebox.showwarning("Warning", "No Data Available")
            else:
                Relaibility_based_proteins.config(text=str((len(NEW_DATA_FRAME["PROTEIN ID"]))))
                #LOC_TREE_var.set(NEW_DATA_FRAME)


        Localization_file_name.set(file_name)
    elif len(single) == 1:
        Localization_combo_box_2.current(0)
    else:
        messagebox.showwarning("Warning","Select Strain")


Localization_combo_box_1_var = StringVar()
Localization_combo_box_1 = AutocompleteCombobox(Localization_Frame, completevalues=short, textvariable=Localization_combo_box_1_var)
Localization_combo_box_1.current(0)
Localization_combo_box_1.bind("<<ComboboxSelected>>", Localization_select)
Localization_combo_box_1.place(relx=0.15, rely=0.45)

Localization_combo_box_Label = Label(Localization_Frame, text="Select Bacteria Genus", bg=Second_Color, fg=Label_Color, font=(font_family, 10))
Localization_combo_box_Label.place(relx=0.15, rely=0)

Localization_combo_box_2_var = StringVar()
Localization_combo_box_2 = ttk.Combobox(Localization_Frame, values=[""], textvariable=Localization_combo_box_2_var)
Localization_combo_box_2.bind("<<ComboboxSelected>>", Localization_select)
Localization_combo_box_2.place(relx=0.3, rely=0.45)

Localization_combo_box_Label_Second = Label(Localization_Frame, text="Select Bacteria Strain", bg=Second_Color, fg=Label_Color, font=(font_family, 10))
Localization_combo_box_Label_Second.place(relx=0.3, rely=0)

Reliability_var = IntVar()
Reliability = ttk.Combobox(Localization_Frame, values=[reliable for reliable in range(10,100, 5)], textvariable=Reliability_var)
Reliability.bind("<<ComboboxSelected>>", Localization_select)
Reliability.place(relx=0.45, rely=0.45)

Reliability_Label = Label(Localization_Frame, text="Reliability Score", bg=Second_Color, fg=Label_Color, font=(font_family, 10))
Reliability_Label.place(relx=0.45, rely=0)

Localization_identity_var = IntVar()
Localization_identity = ttk.Combobox(Localization_Frame, values=[Localization_identity_values for Localization_identity_values in range(10,110, 10)], textvariable=Localization_identity_var)
Localization_identity.bind("<<ComboboxSelected>>", Localization_select)
Localization_identity.place(relx=0.6, rely=0.45)

Localization_identity_Label = Label(Localization_Frame, text="Identity Percentage", bg=Second_Color, fg=Label_Color, font=(font_family, 10))
Localization_identity_Label.place(relx=0.6, rely=0)

Relaibility_based_proteins = Label(Localization_Frame, text="0", bg=Second_Color, font=((8)), fg=Label_Color)
Relaibility_based_proteins.place(relx=0.75, rely=0.1)

Relaibility_based_proteins_Label = Label(Localization_Frame, text="Reliable Proteins", bg=Second_Color, fg=Label_Color)
Relaibility_based_proteins_Label.place(relx=0.75, rely=0.48)

def Localization_Help_btn():
    messagebox.showinfo("Info", "Localization will filter secreted, outer membrane and fimbrium proteins. Select a specific identity percentage. candidates will have equal and greater identity percentage of the selected one.")

Localization_Help = Button(Localization_Frame, image=HELP_ICON, bg=First_Color, borderwidth=0, command=Localization_Help_btn)
Localization_Help.place(relx=0.1, rely=0.15)

Non_Homologs_Frame = Frame(window, width=width_frames, height=60, bg=Second_Color, bd=0)
Non_Homologs_Frame.grid(row=4, sticky="w")

Non_Homologs_Proteins = Label(Non_Homologs_Frame, text="0\nProtein(s)", bg=Second_Color, fg=Label_Color, font=(font_family, 12))
Non_Homologs_Proteins.place(relx=0.87, rely=0.15)

Non_Homologs_var = StringVar()
Non_Homologs_check_btn = Checkbutton(Non_Homologs_Frame, variable=Non_Homologs_var,text="Non Host Homologs" ,onvalue="True", offvalue="False", bg=Second_Color, activebackground=Second_Color, font=(font_family, 11))
Non_Homologs_check_btn.deselect()
Non_Homologs_check_btn.place(relx=0.015, rely=0.4)

Non_Homologs_identity_var = IntVar()
Non_Homologs_identity = ttk.Combobox(Non_Homologs_Frame, values=[Non_Homologs_identity_values for Non_Homologs_identity_values in range(5,55,5)], textvariable=Non_Homologs_identity_var)
Non_Homologs_identity.place(relx=0.2, rely=0.4)

Non_Homologs_identity_Label = Label(Non_Homologs_Frame, text="Identity Percentage", bg=Second_Color, fg=Label_Color, font=(font_family, 10))
Non_Homologs_identity_Label.place(relx=0.2, rely=0)

Non_Homologs_perc_var = IntVar()
Non_Homologs_perc = ttk.Combobox(Non_Homologs_Frame, values=[Non_Homologs_perc_values for Non_Homologs_perc_values in range(5,105,5)], textvariable=Non_Homologs_perc_var)
Non_Homologs_perc.place(relx=0.35, rely=0.4)

Non_Homologs_perc_Label = Label(Non_Homologs_Frame, text="Non Homology Percentage", bg=Second_Color, fg=Label_Color, font=(font_family, 10))
Non_Homologs_perc_Label.place(relx=0.35, rely=0.0)

def Non_Homologs_Help_Btn():
    messagebox.showinfo("Info", "Select specific identity percentage, candidates will have equal and less identity percentage of the selected one.")

Non_Homologs_Help = Button(Non_Homologs_Frame, image=HELP_ICON, bg=First_Color, borderwidth=0, command=Non_Homologs_Help_Btn)
Non_Homologs_Help.place(relx=0.145, rely=0.15)

VFDB_Frame = Frame(window, width=width_frames, height=60, bg=Second_Color, bd=0)
VFDB_Frame.grid(row=5, sticky="w")

PVCs = Label(VFDB_Frame, text="0\nPVC(s)", bg=Second_Color, fg=Label_Color, font=(font_family, 13))
PVCs.place(relx=0.8, rely=0.15)

VFDB_Proteins = Label(VFDB_Frame, text="0\nProtein(s)", bg=Second_Color, fg=Label_Color, font=(font_family, 12))
VFDB_Proteins.place(relx=0.87, rely=0.15)

VFDB_var = StringVar()
VFDB_check_btn = Checkbutton(VFDB_Frame, variable=VFDB_var,text="Virulence Factors",onvalue="True", offvalue="False", bg=Second_Color, activebackground=Second_Color, font=(font_family, 11))
VFDB_check_btn.deselect()
VFDB_check_btn.place(relx=0.015, rely=0.4)

VFDB_identity_var = IntVar()
VFDB_identity = ttk.Combobox(VFDB_Frame, values=[VFDB_identity_values for VFDB_identity_values in range(10,110, 10)], textvariable=VFDB_identity_var)
VFDB_identity.place(relx=0.18, rely=0.4)

VFDB_identity_Label = Label(VFDB_Frame, text="Identity Percentage", bg=Second_Color, fg=Label_Color, font=(font_family, 10))
VFDB_identity_Label.place(relx=0.18, rely=0)

def VFDB_Help_Btn():
    messagebox.showinfo("Info", "Select a specific identity percentatge, Virulence factors will have equal and greater identity percentage of the selected one.")

VFDB_Help = Button(VFDB_Frame, image=HELP_ICON, bg=First_Color, borderwidth=0, command=VFDB_Help_Btn)
VFDB_Help.place(relx=0.133, rely=0.15)

Epitope_Frame = Frame(window, width=width_frames, height=60, bg=Second_Color, bd=0)
Epitope_Frame.grid(row=6, sticky="w")

Epitope_peptides = Label(Epitope_Frame, text="0\nEpitope(s)", bg=Second_Color, fg=Label_Color, font=(font_family, 12))
Epitope_peptides.place(relx=0.87, rely=0.15)

Epitope_var = StringVar()
Epitope_check_btn = Checkbutton(Epitope_Frame, variable=Epitope_var,text="Epitope Mapping" , onvalue="True", offvalue="False", bg=Second_Color, activebackground=Second_Color, font=(font_family, 11))
Epitope_check_btn.deselect()
Epitope_check_btn.place(relx=0.015, rely=0.4)

Epitope_window_B_cell_Label = Label(Epitope_Frame, text="B Cell Length", bg=Second_Color, fg=Label_Color, font=(font_family, 10))
Epitope_window_B_cell_Label.place(relx=0.18, rely=0)

Epitope_window_B_cell_var = StringVar()
Epitope_window_B_cell = ttk.Combobox(Epitope_Frame, width=10, values=[12, 13, 14, 15, 16, 17, 18, 19, 20, "ALL"], textvariable=Epitope_window_B_cell_var)
Epitope_window_B_cell.place(relx=0.18, rely=0.4)

Epitope_identity_Label_B_cell = Label(Epitope_Frame, text="B CELL Identity", bg=Second_Color, fg=Label_Color, font=(font_family, 10))
Epitope_identity_Label_B_cell.place(relx=0.28, rely=0)

Epitope_identity_B_cell_var = IntVar()
Epitope_identity_B_cell = ttk.Combobox(Epitope_Frame, width=10, values=[Epitope_identity_B_cell_values for Epitope_identity_B_cell_values in range(10,110, 10)], textvariable=Epitope_identity_B_cell_var)
Epitope_identity_B_cell.place(relx=0.28, rely=0.4)

Epitope_window_CD8_Label = Label(Epitope_Frame, text="CD8 Length", bg=Second_Color, fg=Label_Color, font=(font_family, 10))
Epitope_window_CD8_Label.place(relx=0.38, rely=0)

Epitope_window_CD8_var = StringVar()
Epitope_window_CD8 = ttk.Combobox(Epitope_Frame, width=10, values=[8, 9, 10, 11, 12, 13, "ALL"], textvariable=Epitope_window_CD8_var)
Epitope_window_CD8.place(relx=0.38, rely=0.4)

Epitope_identity_Label_CD8 = Label(Epitope_Frame, text="CD8 Identity", bg=Second_Color, fg=Label_Color, font=(font_family, 10))
Epitope_identity_Label_CD8.place(relx=0.48, rely=0)

Epitope_identity_CD8_var = IntVar()
Epitope_identity_CD8 = ttk.Combobox(Epitope_Frame, width=10, values=[Epitope_identity_CD8_values for Epitope_identity_CD8_values in range(10,110, 10)], textvariable=Epitope_identity_CD8_var)
Epitope_identity_CD8.place(relx=0.48, rely=0.4)

Epitope_window_CD4_Label = Label(Epitope_Frame, text="CD4 Length", bg=Second_Color, fg=Label_Color, font=(font_family, 10))
Epitope_window_CD4_Label.place(relx=0.58, rely=0)

Epitope_window_CD4_var = StringVar()
Epitope_window_CD4 = ttk.Combobox(Epitope_Frame, width=10, values=[13, 14, 15, 16, 17, 18, 19, 20, "ALL"], textvariable=Epitope_window_CD4_var)
Epitope_window_CD4.place(relx=0.58, rely=0.4)

Epitope_identity_Label_CD4 = Label(Epitope_Frame, text="CD4 Identity", bg=Second_Color, fg=Label_Color, font=(font_family, 10))
Epitope_identity_Label_CD4.place(relx=0.68, rely=0)

Epitope_identity_CD4_var = IntVar()
Epitope_identity_CD4 = ttk.Combobox(Epitope_Frame, width=10, values=[Epitope_identity_CD4_values for Epitope_identity_CD4_values in range(10,110, 10)], textvariable=Epitope_identity_CD4_var)
Epitope_identity_CD4.place(relx=0.68, rely=0.4)

def Epitope_Help_Btn():
    messagebox.showinfo("Info", "Select a specific window length and identity percentage. Epitopes will have equal and greater identity percentage of the selected one.")

Epitope_Help = Button(Epitope_Frame, image=HELP_ICON, bg=First_Color, borderwidth=0, command=Epitope_Help_Btn)
Epitope_Help.place(relx=0.13, rely=0.15)

SUBMIT_BTN_FRAME = Frame(window, width=width_frames, height=50, bg=Second_Color, bd=0)
SUBMIT_BTN_FRAME.grid(row=7, sticky="w")

Fasta_Entry_Btn = Button(SUBMIT_BTN_FRAME, image=SUBMIT_BUTTON,bd=0, command=lambda: threading.Thread(target=FASTA_SUBMIT_BTN).start())
Fasta_Entry_Btn.place(relx=0.4, rely=0)

DEFAULT_VALUES_var = StringVar()
DEFAULT_VALUES = Checkbutton(Must_Evalutate_Frame, variable=DEFAULT_VALUES_var,text="Set Default Values",onvalue="True", offvalue="False", bg=Second_Color, activebackground=Second_Color, font=(font_family, 11), command=SET_DEFAULT_VALUES)
DEFAULT_VALUES.deselect()
DEFAULT_VALUES.place(relx=0.8, rely=0.1)

Footer_Frame = Frame(window, width=width_frames, height=60, bg=Buttons_Colors, bd=0)
Footer_Frame.grid(row=8, sticky="w")

Footer_Frame_Label = Label(Footer_Frame, text="Integrative Biology Laboratory, Atta-ur-Rahman School of Applied Biosciences (ASAB)", bg=Buttons_Colors, fg=Second_Color, font=(font_family, 12))
Footer_Frame_Label.place(relx=0.185, rely=0.05)

Footer_Frame_Label_second = Label(Footer_Frame, text="National University of Sciences and Technology (NUST), Islamabad.", bg=Buttons_Colors, fg=Second_Color, font=(font_family, 12))
Footer_Frame_Label_second.place(relx=0.25, rely=0.43)

Localization_combo_box_1_var.set("All")
Localization_identity.set(60)
Reliability.set(50)
Non_Homologs_identity.set(35)
Non_Homologs_perc.set(50)
VFDB_identity.set(70)
Epitope_identity_CD4.set(40)
Epitope_identity_CD8.set(50)
Epitope_window_CD4.set(20)
Epitope_window_CD8.set(9)
Epitope_window_B_cell.set(15)
Epitope_identity_B_cell.set(50)

window.mainloop()
