import numpy as np
import matplotlib.pyplot as plt
import ROOT
import uproot
import os
from scipy.optimize import curve_fit
import logging
import pandas as pd
import traceback
from scipy.signal import savgol_filter
import pickle

logging.basicConfig(
    filename="../logs/processing_log.txt",
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

def gaussian(x, a, x0, sigma):
    return a * np.exp(-((x - x0) ** 2) / (2 * sigma ** 2))

def calculate_fwhm(amplitude, sigma):
    return 2.355 * abs(sigma)

def analyze_histogram(bin_centers, hist_values, maxfev=10000):
    try:
        hist_values_smooth = savgol_filter(hist_values, 
                                         window_length=min(7, len(hist_values) if len(hist_values) % 2 == 1 else len(hist_values) - 1),
                                         polyorder=3)
        
        peak_idx = np.argmax(hist_values_smooth)
        peak_height = hist_values_smooth[peak_idx]
        peak_pos = bin_centers[peak_idx]
        
        # Estimate initial sigma
        half_height = peak_height / 2
        above_half = hist_values_smooth >= half_height
        if np.sum(above_half) >= 2:
            width_indices = np.where(above_half)[0]
            sigma_estimate = max((bin_centers[width_indices[-1]] - bin_centers[width_indices[0]]) / 2.355, 1e-6)
        else:
            sigma_estimate = max((bin_centers[-1] - bin_centers[0]) / 10, 1e-6)
        
        # Set bounds with minimum separation
        MIN_SEPARATION = 1e-6
        
        lower_bounds = [
            0,  # amplitude cannot be negative
            min(bin_centers),  # mean must be within data range
            sigma_estimate / 100  # sigma must be positive but reasonable
        ]
        
        upper_bounds = [
            peak_height * 2 + MIN_SEPARATION,  # amplitude can be larger than max
            max(bin_centers) + MIN_SEPARATION,  # mean must be within data range
            sigma_estimate * 100 + MIN_SEPARATION  # sigma must be reasonable
        ]
        
        # Ensure bounds are properly separated
        for i in range(len(lower_bounds)):
            if lower_bounds[i] >= upper_bounds[i]:
                lower_bounds[i] = min(lower_bounds[i], upper_bounds[i]) - MIN_SEPARATION
                upper_bounds[i] = max(lower_bounds[i], upper_bounds[i]) + MIN_SEPARATION
        
        p0 = [peak_height, peak_pos, sigma_estimate]
        
        popt, pcov = curve_fit(
            gaussian,
            bin_centers,
            hist_values_smooth,
            p0=p0,
            bounds=(lower_bounds, upper_bounds),
            maxfev=maxfev,
            method='trf'
        )
        
        amplitude, mean, sigma = popt
        fwhm = calculate_fwhm(amplitude, sigma)
        
        # Validate results
        if not (0 < fwhm < (bin_centers[-1] - bin_centers[0])):
            logging.warning("FWHM outside reasonable range, fit may be unreliable")
            return None, None, None
            
        return amplitude, mean, fwhm
        
    except Exception as e:
        logging.error(f"Error in Gaussian fitting: {str(e)}\n{traceback.format_exc()}")
        return None, None, None

def process_root_ttree(root_file, file_name,title = "", output_path = "", summary_data = [], num_bins=100):
    
    if title == "LEAD":
        canvas = root_file.Get("Canvas_1;1")
        if not canvas:
            logging.error(f"TCanvas not found in {file_name}")
            return
        
        primitives = canvas.GetListOfPrimitives()
        histogram = primitives.FindObject("htemp")
        if not histogram:
            logging.error(f"Histogram 'htemp' not found in {file_name}")
            return
        
        # Extract histogram data
        hist_values = [histogram.GetBinContent(i) for i in range(1, histogram.GetNbinsX() + 1)]
        bin_edges = [histogram.GetBinLowEdge(i) for i in range(1, histogram.GetNbinsX() + 2)]
        bin_centers = [(bin_edges[i] + bin_edges[i + 1]) / 2 for i in range(len(bin_edges) - 1)]
        hist_min = bin_edges[0]  # First bin's lower edge
        hist_max = bin_edges[-1]
        entries = histogram.GetEntries()
        if entries > 0:
            mean = histogram.GetMean()
            std_dev = histogram.GetStdDev()
        else:
            mean = std_dev = 0
    else:
        # For TTree
        branch = root_file["mu_14cm;1/Energy__MeV_"]
        energy_data = branch.array(library="np")
        entries = len(energy_data)
        if entries > 0:
            mean = np.mean(energy_data)
            std_dev = np.std(energy_data)
        else:
            mean = std_dev = 0
        hist_min, hist_max = np.percentile(energy_data, [0.1, 99.9])
        bin_edges = np.linspace(hist_min, hist_max, num_bins + 1)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        hist_values, _ = np.histogram(energy_data, bins=bin_edges)
        
    amplitude, _, fwhm = analyze_histogram(bin_centers, hist_values)
    compton_edge = 0.478
    # Save results
    results = {
        "File": file_name,
        "Peak Intensity": amplitude,
        "Mean (Peak Position)": mean,
        "FWHM": fwhm,
        "Entries": entries,
        "Std Dev.": std_dev,
        "Compton Edge": compton_edge
    }
    summary_data.append(results)
    try:

        # Create enhanced plot
        fig = plt.figure(figsize=(12, 8))
        plt.bar(bin_centers, hist_values, width=(bin_edges[1] - bin_edges[0]), 
            color='blue', alpha=0.7, label="Raw Data")
        
        if all(v is not None for v in [amplitude, mean, fwhm]):
            x_fit = np.linspace(hist_min, hist_max, 1000)
            plt.plot(x_fit, gaussian(x_fit, amplitude, mean, fwhm/2.355), 
                    color="red", label="Gaussian Fit", linewidth=2)
        
        if compton_edge is not None:
            plt.axvline(compton_edge, color="green", linestyle="--", 
                    label="Compton Edge")

        plt.title(f"{title} Energy Spectrum - {file_name}", fontsize=14)
        plt.xlim(0.55, 0.75)
        plt.xlabel("Energy (MeV)", fontsize=12)
        plt.ylabel("Counts", fontsize=12)
        plt.legend(fontsize=10)
        plt.grid(True, alpha=0.3)
        
        # Fixed text box formatting
        text_str = [
            f"Peak: {amplitude:.2f}" if isinstance(amplitude, (int, float)) else "Peak: N/A",
            f"FWHM: {fwhm:.4f}" if isinstance(fwhm, (int, float)) else "FWHM: N/A",
            f"Entries: {entries:}",
            f"Mean: {mean:.6f}",
            f"Std Dev: {std_dev:.7f}"
        ]
        text_str = '\n'.join(text_str)
        
        plt.text(1.10,1.10, text_str, transform=plt.gca().transAxes,
                verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
        if output_path != "":    
            save_matplotlib_figure(fig, output_path, file_name)
            plt.close()            
        else:
            plt.show()
        logging.info(f"Successfully processed {file_name}")

    except Exception as e:
        logging.error(f"Error processing file {file_name}: {str(e)}\n{traceback.format_exc()}")

def save_matplotlib_figure(fig, output_path, file_name):
    fig_file = os.path.join(output_path, f"{os.path.splitext(file_name)[0]}.fig")
    with open(fig_file, "wb") as f:
        pickle.dump(fig, f)


def load_fig(fig_file):
    with open(fig_file, "rb") as f:
        fig = pickle.load(f)
    fig.show()

def process_root_file(file_path, file_name, title = "", output_path = "", summary_data = [], num_bins=100):
    """
    Processes a ROOT file to analyze both TTree and TCanvas objects.

    Parameters:
    - file_path: Path to the ROOT file.
    - file_name: Name of the file being processed.
    - output_path: Directory to save output files.
    - summary_data: List to store summary data for both TTree and TCanvas.
    - num_bins: Number of bins for histograms (default: 100).
    """
    if "LEAD" in file_path:
        root_file = ROOT.TFile.Open(file_path)
        title = "LEAD"
    else:
        root_file = uproot.open(file_path)
    if not root_file:
        logging.error(f"Unable to open ROOT file: {file_path}")
        return
    try:
        process_root_ttree(root_file, file_name,title, output_path, summary_data, num_bins)
    except Exception as e:
        logging.error(f"Error processing ROOT file {file_name}: {str(e)}\n{traceback.format_exc()}")
    finally:
        if "LEAD" in file_path:
            root_file.Close()
        else:
            root_file.close()

def main():
    try:
        # Main directory and output paths
        directories = ["../data/LEAD/","../data/ROOT_30MNB5/","../data/ROOT_HY/"]
        for directory_path in directories:
            if not os.path.exists(directory_path):
                raise FileNotFoundError(f"Directory not found: {directory_path}")

            output_path = "../outputs/" + directory_path.removeprefix('../data/').rstrip('/') + "_OUTPUT/"
            os.makedirs(output_path, exist_ok=True)

            # Process all files in the directory
            files = [f for f in os.listdir(directory_path) 
                if os.path.isfile(os.path.join(directory_path, f)) and f.endswith(".root")]

            if not files:
                logging.warning(f"No ROOT files found in {directory_path}")
                return

            summary_data = []

            for file in files:
                full_path = os.path.join(directory_path, file)
                process_root_file(file_path=full_path,file_name=file,output_path=output_path,summary_data=summary_data)
            # Save summary results to a CSV file
            if summary_data:
                summary_df = pd.DataFrame(summary_data)
                summary_csv_path = os.path.join(output_path, "summary_results.csv")
                summary_df.to_csv(summary_csv_path, index=False)
                logging.info(f"Successfully saved summary results to {summary_csv_path}")
            else:
                logging.warning("No summary data was generated")

    except Exception as e:
        logging.error(f"Main execution error: {str(e)}\n{traceback.format_exc()}")

if __name__ == "__main__":
    main()