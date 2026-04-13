import logging
import os
import sys
from datetime import datetime

def setup_logger(module_name: str, log_dir: str = "logs") -> logging.Logger:
    """
    Sets up a professional logging infrastructure that outputs to both
    the console and a timestamped log file.
    
    Args:
        module_name (str): Name of the module/script being executed.
        log_dir (str): Directory where log files will be saved.
        
    Returns:
        logging.Logger: Configured logger object.
    """
    # Create logs directory if it doesn't exist
    os.makedirs(log_dir, exist_ok=True)
    
    # Generate log filename with current date (e.g., pipeline_20231025.log)
    date_str = datetime.now().strftime("%Y%m%d")
    log_file = os.path.join(log_dir, f"pipeline_{date_str}.log")
    
    logger = logging.getLogger(module_name)
    logger.setLevel(logging.INFO)
    
    # Prevent adding handlers multiple times if instantiated repeatedly
    if not logger.handlers:
        # Formatter: [Timestamp] - [Module] - [Level] - Message
        formatter = logging.Formatter(
            '%(asctime)s - [%(name)s] - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        
        # File Handler (Appends to the daily log file)
        file_handler = logging.FileHandler(log_file, mode='a', encoding='utf-8')
        file_handler.setFormatter(formatter)
        file_handler.setLevel(logging.INFO)
        
        # Console Handler (Prints to terminal)
        console_handler = logging.StreamHandler(sys.stdout)
        console_handler.setFormatter(formatter)
        console_handler.setLevel(logging.INFO)
        
        logger.addHandler(file_handler)
        logger.addHandler(console_handler)
        
    return logger
