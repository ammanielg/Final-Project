import os
import numpy as np
import matplotlib.pyplot as plt
import ROOT
import logging

logging.basicConfig(
    filename="../logs/gamma_plot_log.txt",
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S"
)

def extract_gamma_data_from_tcanvas(file_path, canvas_name="Canvas_1"):
    """
    Extract gamma spectrum data from a TCanvas in a ROOT file.

    Parameters:
    - file_path: Path to the ROOT file.
    - canvas_name: Name of the TCanvas object (default: "Canvas_1").

    Returns:
    - bin_centers: Array of bin centers.
    - hist_values: Array of histogram values.
    """
    try:
        root_file = ROOT.TFile.Open(file_path)
        if not root_file or root_file.IsZombie():
            logging.error(f"Unable to open ROOT file: {file_path}")
            return None, None

        # Get the TCanvas object
        canvas = root_file.Get(canvas_name)
        if not canvas:
            logging.error(f"TCanvas '{canvas_name}' not found in the ROOT file.")
            root_file.Close()
            return None, None

        # Iterate over primitives in the canvas to find a histogram
        for primitive in canvas.GetListOfPrimitives():
            if isinstance(primitive, ROOT.TH1):  # Check if the primitive is a histogram
                # Extract histogram data
                bin_edges = np.array([primitive.GetBinLowEdge(i) for i in range(1, primitive.GetNbinsX() + 2)])
                bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
                hist_values = np.array([primitive.GetBinContent(i) for i in range(1, primitive.GetNbinsX() + 1)])

                root_file.Close()
                return bin_centers, hist_values

        logging.error(f"No histogram found in TCanvas '{canvas_name}' for file {file_path}.")
        root_file.Close()
        return None, None

    except Exception as e:
        logging.error(f"Error extracting data from TCanvas in file {file_path}: {e}")
        return None, None


def extract_gamma_data(file_path):
    """
    Extract gamma spectrum data from a ROOT file.

    Parameters:
    - file_path: Path to the ROOT file.

    Returns:
    - bin_centers: Array of bin centers.
    - hist_values: Array of histogram values.
    """
    try:
        branch_name = "Energy__MeV_"
        root_file = ROOT.TFile.Open(file_path)
        if not root_file or root_file.IsZombie():
            logging.error(f"Unable to open ROOT file: {file_path}")
            return None, None

        tree_name = None
        for key in root_file.GetListOfKeys():
            obj = key.ReadObj()
            if isinstance(obj, ROOT.TTree):
                tree_name = obj.GetName()
                break

        if not tree_name:
            logging.error("No TTree found in the ROOT file.")
            return None,None

        tree = root_file.Get(tree_name)
        if not tree:
            logging.error(f"Failed to get TTree: {tree_name}")
            return None,None

        if not tree.GetBranch(branch_name):
            logging.error(f"Branch '{branch_name}' not found in the TTree.")
            return None,None
        # Extract data from the branch
        data = []
        for entry in tree:
            data.append(entry.Energy__MeV_)

        data = np.array(data, dtype=float)
        data = data[~np.isnan(data) & ~np.isinf(data)]  # Clean invalid values

        # Generate histogram
        hist_min, hist_max = np.percentile(data, [0.1, 99.9])
        bin_edges = np.linspace(hist_min, hist_max, 100)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        hist_values, _ = np.histogram(data, bins=bin_edges)

        root_file.Close()
        return bin_centers, hist_values

    except Exception as e:
        logging.error(f"Error extracting data from file {file_path}: {e}")
        return None, None

def plot_gamma_spectrum(directory_path, material, output_path):
    """
    Plot gamma spectra for all thicknesses in a directory.

    Parameters:
    - directory_path: Path to the directory containing ROOT files.
    - material: Name of the material (e.g., 'Lead', 'HY80', '30MnB5').
    - output_path: Path to save the combined plot.
    """
    try:
        files = [f for f in os.listdir(directory_path) if f.endswith(".root")]
        files.sort(key=lambda x: int(x.split("cm")[0]))  # Sort files by thickness

        plt.figure(figsize=(12, 8))

        for file_name in files:
            thickness = file_name.split(".")[0]  # Extract thickness (e.g., '1cm')
            file_path = os.path.join(directory_path, file_name)

            # Extract gamma data
            if material != "Lead":
                bin_centers, hist_values = extract_gamma_data(file_path)
            else:
                bin_centers, hist_values = extract_gamma_data_from_tcanvas(file_path)
            if bin_centers is not None and hist_values is not None:
                plt.plot(bin_centers, hist_values, label=f"{thickness}", alpha=0.7)

        plt.title(f"Gamma Spectrum for {material} at Different Thicknesses")
        plt.xlabel("Energy (MeV)")
        plt.ylabel("Counts")
        plt.legend(title="Thickness")
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.show()
        # Save the combined plot
        output_file = os.path.join(output_path, f"{material}_gamma_spectrum.png")
        plt.savefig(output_file)
        plt.close()

        logging.info(f"Gamma spectrum plot saved: {output_file}")

    except Exception as e:
        logging.error(f"Error plotting gamma spectrum: {e}")

def main():
    materials = {
        "Lead": "../data/LEAD",
        "HY80": "../data/ROOT_HY",
        "30MnB5": "../data/ROOT_30MNB5"
    }
    output_path = "../outputs"
    os.makedirs(output_path, exist_ok=True)

    for material, directory in materials.items():
        plot_gamma_spectrum(directory, material, output_path)

if __name__ == "__main__":
    main()
