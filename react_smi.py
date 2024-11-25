from rdkit import Chem
from rdkit.Chem import AllChem

class Reaction():
    def __init__(self):
        reaction = ""
        deps1 = []
        deps2 = []
        pains = []
        react1 = None
        react2 = None
        products = []
        product_smi = ""

    def set_reaction(self, reaction_text):
        self.reaction = AllChem.ReactionFromSmarts(reaction_text)

    def set_protection_pains(self, deps1=[],deps2=[],pains=[]):
        self.deps1 = [Chem.MolFromSmarts(dep) for dep in deps1]
        self.deps2 = [Chem.MolFromSmarts(dep) for dep in deps2]
        self.pains = [Chem.MolFromSmarts(pain) for pain in pains]

    def filt_product(self,only1=False):
        product_smi = ""
        for product in self.products:
            try:
                mol = product[0]
                product_smi = Chem.MolToSmiles(mol)
                # mol = Chem.MolFromSmiles(Chem.MolToSmiles(product[0]))
                # product_smi = Chem.MolToSmiles(mol)
                count = count + 1
                for pain in self.pains:
                    match = mol.GetSubstructMatches(pain)
                    if match:
                        print(f"fail smiles: {product_smi}")
                        product_smi = ""
                if only1 and count > 1:
                    print(f"fail smiles: {product_smi}")
                    product_smi = ""
            except:
                continue
        return product_smi

    def run_reaction(self, react1="", react2=""):
        if not react1:
            return ""
        elif not react2: # react1 but not react2
            self.react1 = Chem.MolFromSmiles(react1)
            for dep in self.deps1:
                for match in self.react1.GetSubstructMatches(dep):
                    self.react1.GetAtomWithIdx(match[0]).SetProp('_protected','1')
            self.products = self.reaction.RunReactants((self.react1,))
            self.product_smi = self.filt_product(only1=True)
            return self.product_smi
        else: # react1 and react2
            self.react1 = Chem.MolFromSmiles(react1)
            self.react2 = Chem.MolFromSmiles(react2)
            for dep in self.deps1:
                for match in self.react1.GetSubstructMatches(dep):
                    self.react1.GetAtomWithIdx(match[0]).SetProp('_protected','1')
            for dep in self.deps2:
                for match in self.react2.GetSubstructMatches(dep):
                    self.react2.GetAtomWithIdx(match[0]).SetProp('_protected','1')
            self.products = self.reaction.RunReactants((self.react1,self.react2))
            self.product_smi = self.filt_product(only1=True)
            return self.product_smi


if __name__ == "__main__":
    # urea_reaction = Reaction()
    # urea_reaction.set_reaction(reaction_text='[N:1]=[C:2]=[O:3].[N:4]>>[N:1][C:2](=[O:3])[N:4]')
    # urea_reaction.set_protection_pains(
    #     deps1 = [],
    #     deps2 = ['[N;$(N~[N,O,S])]','[N;$(NC~[N,O,S])]'],
    #     pains = ['[C]-[Cl,Br,I]','[C](=[O])-[Cl,Br,I]','[S](=[O])-[Cl,Br,I]']
    # )

    # react1 = "O=C=Nc1cccc(Cl)c1"
    # react2 = "N"
    # products = urea_reaction.run_reaction(react1, react2)
    # print(products)

    tetrazole_reaction = Reaction()
    # tetrazole_reaction.set_reaction(reaction_text='[C:1]#[N:2].[N:3]=[N:4]=[N:5]>>[C:1]1=[N:2][N:3][N:4]=[N:5]1')
    tetrazole_reaction.set_reaction(reaction_text='[C:1]#[N:2]>>[C:1]1=[N:2][N][N]=[N]1')
    tetrazole_reaction.set_protection_pains(
        deps1 = [],
        deps2 = [],
        pains = []
    )

    react1 = "CC#N"
    react2 = ""
    products = tetrazole_reaction.run_reaction(react1, react2)
    print(products)
