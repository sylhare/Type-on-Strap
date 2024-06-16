import tkinter as tk
from tkinter import filedialog, messagebox, ttk
import os
from pathlib import Path
from pkg_resources import resource_filename
import subprocess
import pandas as pd
import threading
import glob
from PIL import ImageTk,Image

class ScrollableFrame(tk.Frame):
    def __init__(self, master, **kwargs):
        tk.Frame.__init__(self, master, **kwargs)
        self.canvas = tk.Canvas(self)
        self.scrollbar = ttk.Scrollbar(self, orient="vertical", command=self.canvas.yview)
        self.scrollable_frame = tk.Frame(self.canvas)

        self.scrollable_frame.bind("<Configure>", lambda e: self.canvas.configure(scrollregion=self.canvas.bbox("all")))

        self.canvas.create_window((0, 0), window=self.scrollable_frame, anchor="nw")
        self.canvas.configure(yscrollcommand=self.scrollbar.set)

        self.canvas.pack(side="left", fill="both", expand=True)
        self.scrollbar.pack(side="right", fill="y")
    
    def center_elements(self):
        for child in self.scrollable_frame.winfo_children():
            child.pack_configure(anchor="center")

class easyROB(tk.Tk):
    def __init__(self):
        super().__init__()

        # Add icon
        path_icon = Path(resource_filename("robert", "report"))
        path_icon = path_icon.joinpath('Robert_icon.ico')
        try: # the icon folder isn't found in some systems
            self.iconbitmap(f'{path_icon}')
        except:
            pass

        # Add logo
        path_logo = Path(resource_filename("robert", "report"))
        path_logo = path_logo.joinpath('Robert_logo1.jpg')
        robert_logo = ImageTk.PhotoImage(Image.open(path_logo))

        # Create a label for the logo with a gray border
        my_logo = tk.Label(self, image=robert_logo, bd=1, relief=tk.SOLID, highlightbackground="gray")
        my_logo.pack(side=tk.TOP, fill=tk.X, padx=10, pady=10)
        my_logo.image = robert_logo
                                            
        #Add tittle
        self.title("easyROB")
        self.geometry("400x800")
        
        # Create a scrollable frame
        self.scrollable_frame = ScrollableFrame(self)
        self.scrollable_frame.pack(fill="both", expand=True)

        # Label and button for CSV file
        self.label = tk.Label(self.scrollable_frame.scrollable_frame, text="Select a CSV file:")
        self.label.pack(pady=10)
        self.file_button = tk.Button(self.scrollable_frame.scrollable_frame, text="Select file", command=self.select_file)
        self.file_button.pack(pady=10)

        # Label for the --y dropdown menu
        self.y_label = tk.Label(self.scrollable_frame.scrollable_frame, text="Select a column for --y (target value to predict):")
        self.y_label.pack(pady=5)

        # Dropdown menu for --y
        self.y_options = tk.StringVar(self.scrollable_frame.scrollable_frame)
        self.y_dropdown = tk.OptionMenu(self.scrollable_frame.scrollable_frame, self.y_options, "")
        self.y_dropdown.pack(pady=5)

        # Dropdown menu for --type
        self.type_label = tk.Label(self.scrollable_frame.scrollable_frame, text="Type of predictions:")
        self.type_label.pack(pady=5)
        self.type_options = tk.StringVar(self.scrollable_frame.scrollable_frame)
        self.type_options.set("Regression")  # Default selection
        self.type_dropdown = tk.OptionMenu(self.scrollable_frame.scrollable_frame, self.type_options, "Regression", "Classification")
        self.type_dropdown.pack(pady=5)

        # Label for the --names dropdown menu
        self.names_label = tk.Label(self.scrollable_frame.scrollable_frame, text="Select a column for --names (names of the data points):")
        self.names_label.pack(pady=5)

        # Dropdown menu for --names
        self.names_options = tk.StringVar(self.scrollable_frame.scrollable_frame)
        self.names_dropdown = tk.OptionMenu(self.scrollable_frame.scrollable_frame, self.names_options, "")
        self.names_dropdown.pack(pady=5)

        # Label for the --csv_test entry field
        self.csv_test_label = tk.Label(self.scrollable_frame.scrollable_frame, text="Select a file for --csv_test (new predictions, optional):")
        self.csv_test_label.pack(pady=5)

        self.csv_test_file_button = tk.Button(self.scrollable_frame.scrollable_frame, text="Select file", command=self.select_csv_test_file)
        self.csv_test_file_button.pack(pady=5)

        # Label for the --start_from_smiles dropdown menu
        self.start_from_smiles_label = tk.Label(self.scrollable_frame.scrollable_frame, text="Start from SMILES (requires AQME and xTB):")
        self.start_from_smiles_label.pack(pady=5)

        # Dropdown menu for --start_from_smiles
        self.start_from_smiles_options = tk.StringVar(self.scrollable_frame.scrollable_frame)
        self.start_from_smiles_options.set("No")  # Default selection
        self.start_from_smiles_dropdown = tk.OptionMenu(self.scrollable_frame.scrollable_frame, self.start_from_smiles_options, "No", "Yes")
        self.start_from_smiles_dropdown.pack(pady=5)

        # Label and Listbox for the --ignore entry field
        self.ignore_label = tk.Label(self.scrollable_frame.scrollable_frame, text="Select columns to ignore:")
        self.ignore_label.pack(pady=5)
        self.ignore_listbox = tk.Listbox(self.scrollable_frame.scrollable_frame, selectmode=tk.MULTIPLE, height=5, width=40)
        self.ignore_listbox.pack(pady=5)

        # Variable to store the path of the selected file
        self.file_path = ""
        self.csv_test_path = ""

        # Make the "Run ROBERT" button larger
        self.run_button = tk.Button(self.scrollable_frame.scrollable_frame, text="Run ROBERT", command=self.run_robert, height=2, width=20)
        self.run_button.pack(pady=10)

        # Progress bar
        self.progress = ttk.Progressbar(self.scrollable_frame.scrollable_frame, orient="horizontal", length=350, mode="indeterminate")
        
        # Load the CSV file to get all column names
        if self.file_path:
            df = pd.read_csv(self.file_path)
            self.all_columns = list(df.columns)

        # Linking change events in drop-down menus
        self.y_options.trace('w', self.update_y_options)
        self.names_options.trace('w', self.update_names_options)
        self.start_from_smiles_options.trace('w', self.update_start_from_smiles_option)

        # Update the geometry of the scrollable_frame after all elements are added
        self.scrollable_frame.update_idletasks()
        self.scrollable_frame.center_elements()

    def update_y_options(self, *args):
        self.update_ignore_listbox()

    def update_names_options(self, *args):
        self.update_ignore_listbox()

    def update_start_from_smiles_option(self, *args):
        self.update_ignore_listbox()

    def update_ignore_listbox(self):
        # Create a set of selected columns in --y and --names
        y_value = self.y_options.get()
        names_value = self.names_options.get()
        selected_columns = set([y_value, names_value])

        # Exclude additional columns if --start_from_smiles is "Yes"
        if self.start_from_smiles_options.get() == "Yes":
            columns_to_exclude = ['smiles', 'charge', 'mult', 'complex_type', 'geom', 'constraints_atoms',
                                'constraints_dist', 'constraints_angle', 'constraints_dihedral']
            selected_columns.update(columns_to_exclude)

        # Update Listbox with column names for --ignore
        self.ignore_listbox.delete(0, 'end')
        for column in self.all_columns:
            if column not in selected_columns:
                self.ignore_listbox.insert('end', column)
                
    def select_file(self):
        self.file_path = filedialog.askopenfilename(filetypes=[("CSV Files", "*.csv")])

        # You can do something with the selected file, for example, display the path in a label
        self.label.config(text=f"Selected file: {self.file_path}")

        # Load the CSV file to get column names for --y and --names
        if self.file_path:
            df = pd.read_csv(self.file_path)
            self.all_columns = list(df.columns)

            # Update dropdown menus with column names
            self.y_options.set("")  # Clear previous selection
            self.y_dropdown['menu'].delete(0, 'end')
            for column in self.all_columns:
                self.y_dropdown['menu'].add_command(label=column, command=tk._setit(self.y_options, column))

            self.names_options.set("")  # Clear previous selection
            self.names_dropdown['menu'].delete(0, 'end')
            for column in self.all_columns:
                self.names_dropdown['menu'].add_command(label=column, command=tk._setit(self.names_options, column))

            # Update Listbox with column names for --ignore
            self.update_ignore_listbox()

            # Set default value for --start_from_smiles
            self.start_from_smiles_options.set("No")

    def select_csv_test_file(self):
        self.csv_test_path = filedialog.askopenfilename(filetypes=[("CSV Files", "*.csv")])
        # You can do something with the selected file, for example, display the path in a label
        self.csv_test_label.config(text=f"Selected file for --csv_test: {self.csv_test_path}")

    def run_robert(self):
        if not self.file_path or not self.y_options.get() or not self.names_options.get():
            messagebox.showwarning("Warning", "Please select a main CSV file, --y, and --names before running ROBERT.")
            return
        
        # Get values for --y, --names, --csv_test, --start_from_smiles, and --ignore from the GUI entries
        y_value = self.y_options.get()

        names_value = self.names_options.get()
        csv_test_value = self.csv_test_path  # Use the path directly
        start_from_smiles_value = self.start_from_smiles_options.get()
        type_value = "clas" if self.type_options.get() == "Classification" else ""
        ignore_selected_indices = self.ignore_listbox.curselection()
        ignore_value = ",".join([self.ignore_listbox.get(index) for index in ignore_selected_indices])

        # Get the directory of the main CSV file
        csv_directory = os.path.dirname(self.file_path)

        # Get only the file name for --csv_name
        csv_name = os.path.basename(self.file_path)

        # Get only the file name for --csv_test
        csv_test_name = os.path.basename(csv_test_value) if csv_test_value else ""

        # Build the command to run robert with the necessary arguments
        command = f'python -m robert --csv_name "{csv_name}" --y "{y_value}" --names "{names_value}"'

        # Add --csv_test if a value is provided
        if csv_test_value:
            command += f' --csv_test "{csv_test_name}"'

         # Add --type if classification is selected
        if type_value:
            command += f' --type "{type_value}"'

        # Add --aqme if "Yes" is selected for --start_from_smiles
        if start_from_smiles_value == "Yes":
            command += ' --aqme'

        # Add --ignore if columns to ignore are selected
        if ignore_value:
            command += f' --ignore "[{ignore_value}]"'

        # Disable the button during execution
        self.run_button.config(state=tk.DISABLED)

        # Show the progress bar
        self.progress.pack(pady=10)
        self.progress.start()

        # Run ROBERT in a separate thread
        thread = threading.Thread(target=self.run_robert_thread, args=(command, csv_directory))
        thread.start()

    def run_robert_thread(self, command, csv_directory):
        try:
            # Run the subprocess and capture the output
            result = subprocess.run(command, shell=True, cwd=csv_directory, check=False, capture_output=True, text=True)

            # Check the return code explicitly
            if result.returncode != 0:
                # Hide the progress bar
                self.progress.stop()
                self.progress.pack_forget()

                # Show an error message if execution fails
                messagebox.showerror("Error", f"Error while running ROBERT:\n{result.stderr}")
            else:
                # Check the last .dat file for the expected line
                last_dat_file = self.find_last_dat_file(csv_directory)
                if last_dat_file:
                    module_name = os.path.splitext(os.path.basename(last_dat_file))[0]
                    module_name = module_name.rsplit("_", 1)[0]  # Remove "_data" part
                    last_line = self.get_last_line_of_file(last_dat_file)
                    if 'Time' in last_line:
                        # Hide the progress bar
                        self.progress.stop()
                        self.progress.pack_forget()

                        # Show a success message after execution
                        messagebox.showinfo("Success", "ROBERT has been executed successfully.")
                    else:
                        # Hide the progress bar
                        self.progress.stop()
                        self.progress.pack_forget()

                        # Show an error message indicating that "Time" is not present in the last line
                        messagebox.showerror("Error", f"Error while running ROBERT. The {module_name} module shows the following error:\n{last_line}")

        except Exception as e:
            # Hide the progress bar
            self.progress.stop()
            self.progress.pack_forget()

            # Show an error message if an exception occurs
            messagebox.showerror("Error", f"Error while running ROBERT:\n{e}")
        finally:
            # Enable the button after execution
            self.run_button.config(state=tk.NORMAL)

    def find_last_dat_file(self, directory):
        # Recursively search for .dat files in the given directory
        dat_files = glob.glob(os.path.join(directory, '**', '*_data.dat'), recursive=True)

        # Return the last modified .dat file
        return max(dat_files, key=os.path.getmtime, default=None)

    def get_last_line_of_file(self, file_path):
        # Read the last line of the given file
        with open(file_path, 'r', encoding='utf-8') as file:
            lines = file.readlines()
            # Remove empty lines at the end of the file
            lines = [line.rstrip('\n') for line in lines if line.strip()]
            if lines:
                # Get the last line after removing empty lines
                return lines[-1]
            else:
                return "File is empty"

if __name__ == "__main__":
    app = easyROB()
    app.mainloop()


