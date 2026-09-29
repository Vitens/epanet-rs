//! Current selection state. Kept separate from `AppState` mainly for testability.

#[derive(Debug, Clone, Default, PartialEq)]
pub struct Selection {
    pub nodes: Vec<String>,
    pub links: Vec<String>,
}

impl Selection {
    pub fn is_empty(&self) -> bool {
        self.nodes.is_empty() && self.links.is_empty()
    }

    pub fn clear(&mut self) {
        self.nodes.clear();
        self.links.clear();
    }

    pub fn single_node(id: impl Into<String>) -> Self {
        Selection {
            nodes: vec![id.into()],
            links: vec![],
        }
    }

    pub fn single_link(id: impl Into<String>) -> Self {
        Selection {
            nodes: vec![],
            links: vec![id.into()],
        }
    }

    pub fn contains_node(&self, id: &str) -> bool {
        self.nodes.iter().any(|n| n == id)
    }

    pub fn contains_link(&self, id: &str) -> bool {
        self.links.iter().any(|l| l == id)
    }

    /// If exactly one node (and no links) is selected, return its id.
    pub fn as_single_node(&self) -> Option<&str> {
        if self.nodes.len() == 1 && self.links.is_empty() {
            Some(&self.nodes[0])
        } else {
            None
        }
    }

    /// If exactly one link (and no nodes) is selected, return its id.
    pub fn as_single_link(&self) -> Option<&str> {
        if self.links.len() == 1 && self.nodes.is_empty() {
            Some(&self.links[0])
        } else {
            None
        }
    }
}
