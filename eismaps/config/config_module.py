import json
import os
from eismaps.utils.logo import print_logo

def custom_input():
    print("Doppler velocity")
    dv_low = input("Enter the lower limit for Doppler velocity (km/s)> ")
    dv_high = input("Enter the upper limit for Doppler velocity (km/s)> ")
    print("Doppler width")
    lw_low = input("Enter the lower limit for Doppler width (km/s)> ")
    lw_high = input("Enter the upper limit for Doppler width (km/s)> ")
    print("Non-thermal velocity")
    nt_low = input("Enter the lower limit for non-thermal velocity (km/s)> ")
    nt_high = input("Enter the upper limit for non-thermal velocity (km/s)> ")

    return {
        "dv_low": dv_low,
        "dv_high": dv_high,
        "lw_low": lw_low,
        "lw_high": lw_high,
        "nt_low": nt_low,
        "nt_high": nt_high
    }

def get_vals():
    cwd = os.getcwd()
    # Initialize my_config with base_path
    my_config = {"base_path": cwd}

    custom = input("Do you want to customise your output? [y,n]> ")

    if custom.lower() == "y":
        custom_options = custom_input()
        my_config.update(custom_options)  # Append custom options to the initial dictionary
    elif custom.lower() == "n":
        print("Default values will be used.")
        # Optionally, define and append default values here
    else:
        print("Invalid input. Please try again.")
        return get_vals()

    return my_config

def create_config(my_config):
    with open("config.json", "w") as f:
        json.dump(my_config, f, indent=2)
    print("Config file created.")

def check_directory():
    current_directory = os.getcwd()

    print('Data anlysis will be performed in your current working directory.')
    print(f'You are currently in: {current_directory}')

    user_input = input("Is this the correct directory? [y,n]> ")
    if user_input.lower() == "n":
        print("Please run the script from the correct directory.")
        exit()
    elif user_input.lower() == "y":
        print("Continuing with the analysis.")

def main():
    check_directory()
    print_logo()
    my_config = get_vals()
    create_config(my_config)

if __name__ == "__main__":
    main()
