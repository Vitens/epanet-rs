//! Parsers and formatters for EPANET time strings (`HH:MM`, `HH:MM:SS`, unit suffixes).

use crate::error::InputError;

/// Parse a time string with optional time unit suffix.
/// Matches EPANET C `hour()` function behavior exactly:
/// - Supports HH:MM or HH:MM:SS clock time format
/// - Supports decimal time with unit suffixes (SECONDS, MINUTES, HOURS, DAYS)
/// - Supports AM/PM suffixes with validation (hours must be < 13)
/// - Defaults to HOURS when no unit is provided for decimal values
/// - Rounds to nearest second (adds 0.5 before truncating)
///
/// Examples:
///   5           = 5 * 3600 sec
///   5 MINUTES   = 5 * 60   sec
///   13:50       = 13*3600 + 50*60 sec
///   1:50 pm     = (12+1)*3600 + 50*60 sec
pub fn parse_time_str(time_str: &str, unit_or_suffix: Option<&str>) -> Result<usize, InputError> {
    let units = unit_or_suffix.unwrap_or("").to_uppercase();

    // Separate clock time into hrs, min, sec (supports HH:MM or HH:MM:SS)
    let time_parts: Vec<&str> = time_str.split(':').collect();

    // Parse up to 3 components (hours, minutes, seconds)
    let mut y = [0.0_f64; 3];
    for (i, part) in time_parts.iter().enumerate().take(3) {
        y[i] = part
            .parse::<f64>()
            .map_err(|_| InputError::new(format!("Invalid time component: {}", part)))?;
    }

    let hours_as_float = match time_parts.len() {
        // Single numeric value with optional unit
        1 => convert_decimal_time_to_hours(y[0], units.as_str())?,

        // Clock time format (HH:MM or HH:MM:SS)
        2..=3 => {
            let combined_hours = y[0] + y[1] / 60.0 + y[2] / 3600.0;
            apply_am_pm_conversion(combined_hours, units.as_str())?
        }

        _ => return Err(InputError::new("Invalid time format".to_string())),
    };

    // Convert hours to seconds with rounding (matches C: (long)(3600.0 * y + 0.5))
    let seconds = (3600.0 * hours_as_float + 0.5) as usize;

    Ok(seconds)
}

/// Convert a decimal time value to hours based on the unit
fn convert_decimal_time_to_hours(value: f64, unit: &str) -> Result<f64, InputError> {
    match unit {
        "" => Ok(value),
        "SECOND" | "SECONDS" | "SEC" => Ok(value / 3600.0),
        "MINUTE" | "MINUTES" | "MIN" => Ok(value / 60.0),
        "HOUR" | "HOURS" => Ok(value),
        "DAY" | "DAYS" => Ok(value * 24.0),
        "AM" | "PM" => apply_am_pm_conversion(value, unit),
        _ => Err(InputError::new(format!("Invalid time unit: {}", unit))),
    }
}

/// Apply AM/PM conversion to hours value
/// AM: 12.xx → 0.xx (midnight), others stay as-is
/// PM: 12.xx stays as-is (noon), others add 12
fn apply_am_pm_conversion(hours: f64, suffix: &str) -> Result<f64, InputError> {
    match suffix {
        "" => Ok(hours),
        "AM" => {
            if hours >= 13.0 {
                return Err(InputError::new(format!("Invalid hour for AM: {}", hours)));
            }
            Ok(if hours >= 12.0 { hours - 12.0 } else { hours })
        }
        "PM" => {
            if hours >= 13.0 {
                return Err(InputError::new(format!("Invalid hour for PM: {}", hours)));
            }
            Ok(if hours >= 12.0 { hours } else { hours + 12.0 })
        }
        _ => Err(InputError::new(format!("Invalid time suffix: {}", suffix))),
    }
}

pub fn seconds_to_hhmm(seconds: usize) -> String {
    let hours = seconds / 3600;
    let minutes = (seconds % 3600) / 60;
    format!("{:02}:{:02}", hours, minutes)
}

pub fn seconds_to_hhmmss(seconds: usize) -> String {
    let hours = seconds / 3600;
    let minutes = (seconds % 3600) / 60;
    let seconds = seconds % 60;
    format!("{:02}:{:02}:{:02}", hours, minutes, seconds)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hhmm_am_pm() {
        assert_eq!(parse_time_str("12:00", Some("AM")).unwrap(), 0);
        assert_eq!(parse_time_str("12:00", Some("PM")).unwrap(), 12 * 3600);
        assert_eq!(parse_time_str("1:30", Some("AM")).unwrap(), 3600 + 30 * 60);
        assert_eq!(
            parse_time_str("1:30", Some("PM")).unwrap(),
            13 * 3600 + 30 * 60
        );
        assert_eq!(
            parse_time_str("11:59", Some("PM")).unwrap(),
            23 * 3600 + 59 * 60
        );
    }

    #[test]
    fn test_hhmm_no_suffix() {
        assert_eq!(parse_time_str("10:15", None).unwrap(), 10 * 3600 + 15 * 60);
        assert_eq!(parse_time_str("0:45", None).unwrap(), 45 * 60);
    }

    #[test]
    fn test_hhmmss_format() {
        assert_eq!(
            parse_time_str("1:30:45", None).unwrap(),
            (3600.0 * (1.0 + 30.0 / 60.0 + 45.0 / 3600.0) + 0.5) as usize
        );
        assert_eq!(parse_time_str("0:0:30", None).unwrap(), 30);
    }

    #[test]
    fn test_numeric_with_units() {
        assert_eq!(parse_time_str("2", Some("HOURS")).unwrap(), 2 * 3600);
        assert_eq!(parse_time_str("15", Some("MIN")).unwrap(), 15 * 60);
        assert_eq!(parse_time_str("30", Some("SEC")).unwrap(), 30);
        assert_eq!(parse_time_str("1", Some("DAY")).unwrap(), 86_400);
    }

    #[test]
    fn test_numeric_am_pm() {
        assert_eq!(parse_time_str("12", Some("AM")).unwrap(), 0);
        assert_eq!(parse_time_str("1", Some("AM")).unwrap(), 3600);
        assert_eq!(parse_time_str("11", Some("AM")).unwrap(), 11 * 3600);

        assert_eq!(parse_time_str("12", Some("PM")).unwrap(), 12 * 3600);
        assert_eq!(parse_time_str("1", Some("PM")).unwrap(), 13 * 3600);
        assert_eq!(parse_time_str("11", Some("PM")).unwrap(), 23 * 3600);

        assert!(parse_time_str("13", Some("AM")).is_err());
        assert!(parse_time_str("13", Some("PM")).is_err());
    }

    #[test]
    fn test_default_unit_is_hours() {
        assert_eq!(parse_time_str("3", None).unwrap(), 3 * 3600);
    }

    #[test]
    fn test_rounding() {
        assert_eq!(parse_time_str("5.5", None).unwrap(), 19800);
        assert_eq!(parse_time_str("1.999", None).unwrap(), 7196);
        assert_eq!(parse_time_str("2.0003", None).unwrap(), 7201);
    }

    #[test]
    fn test_epanet_c_examples() {
        assert_eq!(parse_time_str("5", None).unwrap(), 5 * 3600);
        assert_eq!(parse_time_str("5", Some("MINUTES")).unwrap(), 5 * 60);
        assert_eq!(parse_time_str("13:50", None).unwrap(), 13 * 3600 + 50 * 60);
        assert_eq!(
            parse_time_str("1:50", Some("PM")).unwrap(),
            (12 + 1) * 3600 + 50 * 60
        );
    }

    #[test]
    fn test_invalid_inputs() {
        assert!(parse_time_str("abc", None).is_err());
        assert!(parse_time_str("10", Some("WEEKS")).is_err());
        assert!(parse_time_str("xx:yy", None).is_err());
    }

    #[test]
    fn test_zero_value() {
        assert_eq!(parse_time_str("0", None).unwrap(), 0);
        assert_eq!(parse_time_str("0", Some("HOURS")).unwrap(), 0);
        assert_eq!(parse_time_str("0", Some("MINUTES")).unwrap(), 0);
        assert_eq!(parse_time_str("0", Some("SECONDS")).unwrap(), 0);
        assert_eq!(parse_time_str("0:00", None).unwrap(), 0);
        assert_eq!(parse_time_str("0:0:0", None).unwrap(), 0);
    }
}
