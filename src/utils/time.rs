//! Parsers and formatters for EPANET time strings (`HH:MM`, `HH:MM:SS`, unit suffixes).

use crate::error::InputError;

/// Parse a time string with optional time unit suffix.
/// Supports formats:
/// - "HH:MM" with optional AM/PM suffix
/// - Numeric value with unit (HOUR, MINUTE, MIN, SECOND, SEC, DAY, AM, PM)
/// - If no unit is provided, defaults to HOURS
pub fn parse_time_str(time_str: &str, unit_or_suffix: Option<&str>) -> Result<usize, InputError> {
    let seconds: usize;

    // If time contains ":", parse as HH:MM format
    if time_str.contains(':') {
        let mut time_parts = time_str.split(':');

        let hours = time_parts
            .next()
            .ok_or_else(|| InputError::new("Missing hours in time"))?
            .parse::<usize>()
            .map_err(|_| InputError::new(format!("Invalid hours value in time: {}", time_str)))?;

        let minutes = time_parts
            .next()
            .ok_or_else(|| InputError::new("Missing minutes in time"))?
            .parse::<usize>()
            .map_err(|_| InputError::new(format!("Invalid minutes value in time: {}", time_str)))?;

        let mut hour24 = hours;

        if let Some(suffix) = unit_or_suffix {
            match suffix.to_uppercase().as_str() {
                "AM" => hour24 = hours % 12,          // 12 AM → 0
                "PM" => hour24 = (hours % 12) + 12,   // 12 PM → 12
                _ => {}
            }
        }

        seconds = hour24 * 3600 + minutes * 60;
    } else {
        // Parse as numeric value (integer or float) with optional time unit
        let value = time_str
            .parse::<f64>()
            .map_err(|_| InputError::new(format!("Invalid time value: {}", time_str)))?;

        let mut unit = unit_or_suffix.unwrap_or("HOURS").to_uppercase();

        // Remove trailing "S" for singular form
        if unit.ends_with('S') && unit != "HOURS" {
            unit.pop();
        }

        seconds = match unit.as_str() {
            "HOUR" | "HOURS" => (value * 3600.0) as usize,
            "MINUTE" | "MIN" => (value * 60.0) as usize,
            "SECOND" | "SEC" => value as usize,
            "DAY" => (value * 86_400.0) as usize,
            "AM" => ((value as usize % 12)) * 3600,
            "PM" => (((value as usize) % 12) + 12) * 3600,
            _ => return Err(InputError::new(format!("Invalid time unit: {}", unit))),
        };
    }

    Ok(seconds)
}

pub fn seconds_to_hhmm(seconds: usize) -> String {
    let hours = seconds / 3600;
    let minutes = (seconds % 3600) / 60;
    format!("{:02}:{:02}", hours, minutes)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_hhmm_am_pm() {
        assert_eq!(parse_time_str("12:00", Some("AM")).unwrap(), 0);
        assert_eq!(parse_time_str("12:00", Some("PM")).unwrap(), 12 * 3600);
        assert_eq!(parse_time_str("1:30", Some("AM")).unwrap(), 1 * 3600 + 30 * 60);
        assert_eq!(parse_time_str("1:30", Some("PM")).unwrap(), 13 * 3600 + 30 * 60);
        assert_eq!(parse_time_str("11:59", Some("PM")).unwrap(), 23 * 3600 + 59 * 60);
    }

    #[test]
    fn test_hhmm_no_suffix() {
        assert_eq!(parse_time_str("10:15", None).unwrap(), 10 * 3600 + 15 * 60);
        assert_eq!(parse_time_str("0:45", None).unwrap(), 45 * 60);
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
        assert_eq!(parse_time_str("12", Some("PM")).unwrap(), 12 * 3600);
        assert_eq!(parse_time_str("1", Some("AM")).unwrap(), 1 * 3600);
        assert_eq!(parse_time_str("1", Some("PM")).unwrap(), 13 * 3600);
    }

    #[test]
    fn test_default_unit_is_hours() {
        assert_eq!(parse_time_str("3", None).unwrap(), 3 * 3600);
    }

    #[test]
    fn test_invalid_inputs() {
        assert!(parse_time_str("abc", None).is_err());
        assert!(parse_time_str("10", Some("WEEKS")).is_err());
        assert!(parse_time_str("xx:yy", None).is_err());
    }
}
