class StatModule:
    def __init__(self, tokenised_genes, strictest_pval, module_id=None):
        self.module_id = module_id
        self.tokenised_genes = tokenised_genes
        self.strictest_pval = strictest_pval
        self.occurences = None
        self.n_genes = len(tokenised_genes)
        self.n_domains = sum(len(gene) for gene in tokenised_genes)

    def __repr__(self):
        return f"StatModule({self.tokenised_genes}, {self.strictest_pval})"

    def __eq__(self, other):
        if not isinstance(other, StatModule):
            return False
        return (self.tokenised_genes == other.tokenised_genes and
                self.strictest_pval == other.strictest_pval)

    def __hash__(self):
        # Convert tokenised_genes to a tuple of tuples to make it hashable
        tokenised_genes_hashable = tuple(tuple(gene) for gene in self.tokenised_genes)
        return hash((tokenised_genes_hashable, self.strictest_pval))

    def to_dict(self):
        """
        Converts the StatModule object to a dictionary representation.

        Returns:
            dict: A dictionary representation of the StatModule object.
        """
        return {
            "module_id": self.module_id,
            "n_genes": self.n_genes,
            "n_domains": self.n_domains,
            "strictest_pval": self.strictest_pval,
            "occurences": self.occurences,
            "tokenised_genes": self.tokenised_genes,
        }
